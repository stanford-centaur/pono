#include "frontends/coreir_encoder.h"
#include "utils/logger.h"

#include <iostream>
#include <set>
#include <sstream>

using namespace CoreIR;
using namespace smt;
using namespace std;

// helpers

bool instance_of(CoreIR::Instance * inst, std::string ns, std::string name)
{
  auto mod = inst->getModuleRef();
  string modns = mod->getNamespace()->getName();
  if (modns == "corebit") {
    // corebit is not generated
    return (ns == "corebit") && (name == mod->getName());
  }
  if (!mod->isGenerated()) {
    return false;
  }
  auto gen = mod->getGenerator();
  return gen->getName() == name && gen->getNamespace()->getName() == ns;
}

// operator map

const unordered_map<string, PrimOp> boolopmap(
    { { "and", And }, { "or", Or }, { "xor", Xor }, { "not", Not } });

const unordered_map<string, PrimOp> bvopmap(
    { { "not", BVNot },   { "and", BVAnd },   { "or", BVOr },
      { "xor", BVXor },   { "shl", BVShl },   { "lshr", BVLshr },
      { "ashr", BVAshr }, { "neg", BVNeg },   { "add", BVAdd },
      { "sub", BVSub },   { "mul", BVMul },   { "udiv", BVUdiv },
      { "urem", BVUrem }, { "sdiv", BVSdiv }, { "srem", BVSrem },
      { "smod", BVSmod }, { "eq", Equal },    { "slt", BVSlt },
      { "sgt", BVSgt },   { "sle", BVSle },   { "sge", BVSge },
      { "ult", BVUlt },   { "ugt", BVUgt },   { "ule", BVUle },
      { "uge", BVUge } });

namespace pono {

// static functions

Module * CoreIREncoder::read_coreir_file(Context * c, std::string filename)
{
  Module * m;
  // load file
  if (!loadFromFile(c, filename, &m)) {
    c->die();
    throw PonoException("Error reading CoreIR file: " + filename);
  }
  return m;
}

// member functions

void CoreIREncoder::encode()
{
  // expecting top_ to be non-null
  assert(top_);

  /* running passes prints to stdout -- redirect output */
  // save old stdout
  streambuf * stdout = cout.rdbuf();
  stringstream ss;
  // redirect
  cout.rdbuf(ss.rdbuf());
  c_->runPasses(passes_, { "global" });
  // replace stdout
  cout.rdbuf(stdout);

  // start processing module
  def_ = top_->getDef();

  // used to determine which inputs of an instance have been processed
  unordered_map<Instance *, set<Wireable *>> covered_inputs;

  // create state elements and store the number of inputs for each instance
  vector<Instance *> instances;
  set<Instance *> state_elements;
  unordered_map<Instance *, size_t> num_inputs;
  bool arst;
  for (auto ipair : def_->getInstances()) {
    type_ = ipair.second->getType();
    if (instance_of(ipair.second, "coreir", "reg")
        || (arst = instance_of(ipair.second, "coreir", "reg_arst"))) {
      if (!ipair.second->getModArgs().at("clk_posedge")->get<bool>()) {
        // cannot abstract clock if there is negative edge clock behavior
        can_abstract_clock_ = false;
      }

      if (arst) {
        // cannot abstract clock if there are asynchronous resets
        can_abstract_clock_ = false;
      }

      state_elements.insert(ipair.second);
      // put registers into instances first (processed last)
      instances.push_back(ipair.second);
    } else if (instance_of(ipair.second, "coreir", "mem")) {
      throw PonoException(
          "Have not implemented memory encoding yet. Please post test case to "
          "pono GitHub issues.");
    }

    size_t n = 0;
    for (auto elem : ipair.second->getSelects()) {
      Type * t = elem.second->getType();
      if (t->isInput() || t->isInOut()) {
        n++;
      }
    }
    num_inputs[ipair.second] = n;
    if (!n) {
      // nodes with no inputs should be processed first
      instances.push_back(ipair.second);
    }
  }

  // create symbols for interface inputs and states for clocks
  for (auto elem : def_->getInterface()->getSelects()) {
    // flip the type (want to view from inside the module)
    type_ = elem.second->getType()->getFlipped();

    // both clock and arst need to be states because we need to call
    // next on them for edge behavior
    bool is_clk;
    if ((is_clk = elem.second->getType()->toString() == "coreir.clk")
        || elem.second->getType()->toString() == "coreir.arst") {
      t_ = ts_.make_statevar(elem.first, boolsort_);
      w2term_[elem.second] = t_;
      if (is_clk) {
        num_clocks_++;
      }
    } else if (type_->isInput() || type_->isInOut()) {
      sort_ = compute_sort(elem.second);
      t_ = ts_.make_inputvar(elem.first, sort_);
      w2term_[elem.second] = t_;
    } else {
      // if not a clock or global input, don't need to connect it yet
      continue;
    }

    // assign all the connections to these global inputs
    Wireable * dst;
    Wireable * parent;
    Instance * parent_inst;
    for (auto conn : elem.second->getLocalConnections()) {
      wire_connection(conn);
      dst = conn.second;
      Type * typ = dst->getType();

      // expecting to have a destination with type input or InOut
      assert(typ->isInput() || typ->isInOut());

      // parent is either an instance or a top-level input
      parent = dst->getTopParent();
      if (Instance::classof(parent)) {
        parent_inst = dyn_cast<Instance>(parent);

        Wireable * dst_parent = dst;
        if (isa<CoreIR::Select>(dst)
            && isNumber(cast<CoreIR::Select>(dst)->getSelStr())) {
          // shouldn't count bit-selects as individual inputs
          // need to get parent
          dst_parent = cast<CoreIR::Select>(dst)->getParent();
        }
        covered_inputs[parent_inst].insert(dst_parent);
        // if all inputs are driven, then add onto stack to be processed
        if (num_inputs[parent_inst] == covered_inputs.at(parent_inst).size()) {
          instances.push_back(parent_inst);
        }
      }
    }
  }

  // can't abstract the clock if there's more than one
  can_abstract_clock_ &= (num_clocks_ <= 1);

  logger.log(1,
             "INFO {} abstract clock for CoreIR Module {}",
             can_abstract_clock_ ? "can" : "cannot",
             top_->getName());

  if (!can_abstract_clock_ && force_abstract_clock_) {
    logger.log(1, "INFO forcing clock abstraction based on user-provided flag");
  }

  // process the rest in topological order
  size_t processed_instances = 0;
  unordered_set<Instance *> visited_instances;
  while (instances.size()) {
    inst_ = instances.back();
    instances.pop_back();

    if (visited_instances.find(inst_) != visited_instances.end()) {
      // this module has already been processed
      continue;
    } else {
      // mark this instance as visited
      visited_instances.insert(inst_);
    }

    Wireable * inst_out = process_instance(inst_);
    processed_instances++;

    // check everything connected to outputs
    // append to list if all other inputs have been seen already

    Wireable * dst;
    Wireable * parent;
    Instance * parent_inst;
    for (Connection conn : inst_out->getLocalConnections()) {
      wire_connection(conn);
      dst = conn.second;
      type_ = dst->getType();

      // expecting to have a destination with type input or InOut
      assert(type_->isInput() || type_->isInOut());

      // parent is either an instance or a top-level input
      parent = dst->getTopParent();
      if (Instance::classof(parent)) {
        parent_inst = dyn_cast<Instance>(parent);

        // state_elements have already been added to instances so ignore those
        if (instance_of(parent_inst, "coreir", "reg")
            || instance_of(parent_inst, "coreir", "reg_arst")) {
          continue;
        }

        Wireable * dst_parent = dst;
        if (isa<CoreIR::Select>(dst)
            && isNumber(cast<CoreIR::Select>(dst)->getSelStr())) {
          // shouldn't count bit-selects as individual inputs
          // need to get parent
          dst_parent = cast<CoreIR::Select>(dst)->getParent();
        }
        covered_inputs[parent_inst].insert(dst_parent);

        // if all inputs are driven, then add onto stack to be processed
        if (num_inputs[parent_inst] == covered_inputs.at(parent_inst).size()
            // optimization -- don't even add to stack if already visited
            && visited_instances.find(parent_inst) == visited_instances.end()) {
          instances.push_back(parent_inst);
        }
      }
    }
  }

  if (processed_instances != def_->getInstances().size()) {
    throw PonoException("Issue: not all instances processed in CoreIR Encoder");
  }

  // now make a second pass over state_elements to assign the next state updates
  for (auto st : state_elements) {
    process_state_element(st);
  }
}

Wireable * CoreIREncoder::process_instance(CoreIR::Instance * inst)
{
  mod_ = inst->getModuleRef();
  string nsname = mod_->getNamespace()->getName();
  string name = mod_->getName();
  t_ = nullptr;
  if (nsname == "corebit" && boolopmap.find(name) != boolopmap.end()) {
    if (name != "not") {
      t_ = solver_->make_term(boolopmap.at(name),
                              w2term_.at(inst->sel("in0")),
                              w2term_.at(inst->sel("in1")));
    } else {
      // special case for not because unary
      t_ = solver_->make_term(boolopmap.at(name), w2term_.at(inst->sel("in")));
    }
  } else if (nsname == "coreir" && bvopmap.find(name) != bvopmap.end()) {
    if (name != "not" && name != "neg") {
      t_ = solver_->make_term(bvopmap.at(name),
                              w2term_.at(inst->sel("in0")),
                              w2term_.at(inst->sel("in1")));
    } else {
      // special case for not and neg because unary
      t_ = solver_->make_term(bvopmap.at(name), w2term_.at(inst->sel("in")));
    }
  } else if (name == "reg" || name == "reg_arst") {
    // NOTE: inputs to state_elements are not wired up until later
    sort_ = solver_->make_sort(
        BV, inst->getModuleRef()->getGenArgs().at("width")->get<int>());
    t_ = ts_.make_statevar(inst->toString(), sort_);
  } else if (nsname == "coreir" && name == "const") {
    size_t w = mod_->getGenArgs().at("width")->get<int>();
    sort_ = solver_->make_sort(BV, w);
    t_ = solver_->make_term(
        (inst->getModArgs().at("value"))->get<BitVec>().binary_string(),
        sort_,
        2);
  } else if (nsname == "corebit" && name == "const") {
    sort_ = boolsort_;
    t_ = solver_->make_term((inst->getModArgs().at("value"))->get<bool>());
  } else if (name == "mux") {
    Term cond = w2term_.at(inst->sel("sel"));
    Term in0 = w2term_.at(inst->sel("in0"));
    Term in1 = w2term_.at(inst->sel("in1"));
    // in1 and in0 swap because a mux selects the first value when sel is 0
    // (e.g. false)
    t_ = solver_->make_term(Ite, cond, in1, in0);
  } else if (name == "slice") {
    Values values = mod_->getGenArgs();
    int hi = values.at("hi")->get<int>();
    int lo = values.at("lo")->get<int>();
    t_ = solver_->make_term(Op(Extract, hi, lo), w2term_.at(inst->sel("in")));
  } else if (name == "concat") {
    t_ = solver_->make_term(
        Concat, w2term_.at(inst->sel("in0")), w2term_.at(inst->sel("in1")));
  } else if (nsname == "coreir" && name == "undriven") {
    sort_ = solver_->make_sort(BV, mod_->getGenArgs().at("width")->get<int>());
    t_ = ts_.make_inputvar(inst->toString(), sort_);
  } else if (nsname == "corebit" && name == "undriven") {
    t_ = ts_.make_inputvar(inst->toString(), boolsort_);
  } else if (name == "andr") {
    // reduce and over bits is only 1 if all bits are 1
    Term in = w2term_.at(inst->sel("in"));
    Sort insort = in->get_sort();
    Term allones =
        solver_->make_term(std::string(insort->get_width(), '1'), insort, 2);
    t_ = solver_->make_term(Equal, in, allones);
  } else if (name == "orr") {
    // reduce or over bits is 1 unless all bits are zero
    Term in = w2term_.at(inst->sel("in"));
    Sort insort = in->get_sort();
    Term zero = solver_->make_term(0, insort);
    t_ = solver_->make_term(Distinct, in, zero);
  } else if (name == "xorr") {
    Term in = w2term_.at(inst->sel("in"));
    Sort insort = in->get_sort();
    int idx = insort->get_width() - 1;
    t_ = solver_->make_term(Op(Extract, idx, idx), in);
    while (idx > 0) {
      idx--;
      t_ = solver_->make_term(
          BVXor, t_, solver_->make_term(Op(Extract, idx, idx), in));
    }
  } else {
    cout << "got instance " << inst->toString() << " : "
         << inst->getModuleRef()->getName() << " but don't know what to do yet!"
         << endl;
  }

  // some modules don't have an output
  if (mod_->getName() == "term") {
    return inst;
  }

  if (!t_) {
    throw PonoException("CoreIREncoder error: no term created for module type: "
                        + mod_->getName());
  }

  w2term_[inst] = t_;
  w2term_[inst->sel("out")] = t_;
  ts_.name_term(inst->sel("out")->toString(), t_);
  return inst->sel("out");
}

void CoreIREncoder::process_state_element(Instance * st)
{
  // expecting a state element
  assert(instance_of(st, "coreir", "reg")
         || instance_of(st, "coreir", "reg_arst")
         || instance_of(st, "coreir", "mem"));

  if (can_abstract_clock_ || force_abstract_clock_) {
    if (w2term_.find(st->sel("in")) != w2term_.end()) {
      ts_.assign_next(w2term_.at(st), w2term_.at(st->sel("in")));
    } else {
      logger.log(1, "Warning: no driver for register {}", st->toString());
    }

    Values vals = st->getModArgs();
    if (vals.find("init") != vals.end()) {
      Term stterm = w2term_.at(st);
      Term initval =
          solver_->make_term(vals.at("init")->get<BitVec>().binary_string(),
                             stterm->get_sort(),
                             2);
      ts_.constrain_init(solver_->make_term(Equal, stterm, initval));
    }

  } else {
    Term clk;
    if (w2term_.find(st->sel("clk")) != w2term_.end()) {
      clk = w2term_.at(st->sel("clk"));
    } else {
      throw PonoException("Clock not wired for register: " + st->toString());
    }

    // expecting clk to be a state variable
    assert(ts_.statevars().find(clk) != ts_.statevars().end());

    Term in;
    if (w2term_.find(st->sel("in")) != w2term_.end()) {
      in = w2term_.at(st->sel("in"));
    } else {
      logger.log(1, "Warning: no driver for register {}", st->toString());
      in = ts_.make_inputvar(st->sel("in")->toString(),
                             compute_sort(st->sel("in")));
    }

    assert(w2term_.find(st) != w2term_.end());
    Term cur_stterm = w2term_.at(st);

    bool clk_posedge = st->getModArgs().at("clk_posedge")->get<bool>();
    Term active_clk;
    if (clk_posedge) {
      active_clk =
          solver_->make_term(And, solver_->make_term(Not, clk), ts_.next(clk));
    } else {
      active_clk =
          solver_->make_term(And, clk, solver_->make_term(Not, ts_.next(clk)));
    }

    Term next_st = solver_->make_term(Ite, active_clk, in, cur_stterm);

    Term initval;
    Values vals = st->getModArgs();
    if (vals.find("init") != vals.end()) {
      initval =
          solver_->make_term(vals.at("init")->get<BitVec>().binary_string(),
                             cur_stterm->get_sort(),
                             2);
    }

    bool async = instance_of(st, "coreir", "reg_arst");
    if (async) {
      bool arst_posedge = st->getModArgs().at("arst_posedge")->get<bool>();
      Term arst_driver;
      if (w2term_.find(st->sel("arst")) != w2term_.end()) {
        arst_driver = w2term_.at(st->sel("arst"));
        if (!ts_.only_curr(arst_driver)) {
          throw PonoException(
              "Driver for ARST has non-state variables -- causes semantic "
              "issues with transition system.");
        }
      } else {
        logger.log(
            1, "Warning: no driver for register arst: {}", st->toString());
        arst_driver = ts_.make_statevar(st->sel("arst")->toString(), boolsort_);
      }

      Term active_arst;
      if (arst_posedge) {
        active_arst = solver_->make_term(
            And, solver_->make_term(Not, arst_driver), ts_.next(arst_driver));
      } else {
        active_arst = solver_->make_term(
            And, arst_driver, solver_->make_term(Not, ts_.next(arst_driver)));
      }

      next_st = solver_->make_term(Ite, active_arst, initval, next_st);
    }

    ts_.constrain_trans(
        solver_->make_term(Equal, ts_.next(cur_stterm), next_st));

    if (initval) {
      ts_.constrain_init(solver_->make_term(Equal, cur_stterm, initval));
    }
  }
}

void CoreIREncoder::wire_connection(Connection conn)
{
  Wireable * src = conn.first;
  Wireable * dst = conn.second;
  // expecting destination to be an input
  // this method should only be called on connections from an output to a
  // destination
  assert(dst->getType()->isInput() || dst->getType()->isInOut());

  bool src_bit_select = isa<CoreIR::Select>(src)
                        && isNumber(cast<CoreIR::Select>(src)->getSelStr());
  bool dst_bit_select = isa<CoreIR::Select>(dst)
                        && isNumber(cast<CoreIR::Select>(dst)->getSelStr());

  Term tmpterm;
  CoreIR::Select * src_sel;
  CoreIR::Select * dst_sel;
  if (src_bit_select && !dst_bit_select) {
    src_sel = cast<CoreIR::Select>(src);
    size_t idx = stoi(src_sel->getSelStr());
    tmpterm = solver_->make_term(Op(Extract, idx, idx), t_);
    tmpterm = solver_->make_term(Equal, tmpterm, bv1_);
    w2term_[src] = tmpterm;
  } else if (!src_bit_select && dst_bit_select) {
    dst_sel = cast<CoreIR::Select>(dst);
    Wireable * parent = dst_sel->getParent();
    size_t idx = stoi(dst_sel->getSelStr());

    Term tparent;
    if (w2term_.find(parent) == w2term_.end()) {
      // create new "input" (actually more of a definition) for dst parent
      // need a forward reference for it
      sort_ = compute_sort(parent);
      tparent = ts_.make_inputvar(parent->toString(), sort_);
      // cache this symbol
      w2term_[parent] = tparent;
    } else {
      tparent = w2term_.at(parent);
    }

    // expecting a bit-vector, cannot select from a Bool
    assert(tparent->get_sort()->get_sort_kind() == BV);

    tmpterm = solver_->make_term(Op(Extract, idx, idx), tparent);
    // would normally expect a boolean
    // but some solvers (e.g. boolector)
    // alias Bool and BV[1], thus we would get BV[1] here
    if (t_->get_sort()->get_sort_kind() == BOOL) {
      // convert to Bool
      tmpterm = solver_->make_term(Equal, tmpterm, bv1_);
    }
    // constrain to be equivalent
    ts_.add_constraint(solver_->make_term(Equal, t_, tmpterm));

  } else if (src_bit_select && dst_bit_select) {
    src_sel = cast<CoreIR::Select>(src);
    dst_sel = cast<CoreIR::Select>(dst);
    Wireable * dst_parent = dst_sel->getParent();
    size_t src_idx = stoi(src_sel->getSelStr());
    size_t dst_idx = stoi(dst_sel->getSelStr());

    Term term_dst_parent;
    if (w2term_.find(dst_parent) == w2term_.end()) {
      // create new "input" (actually more of a definition) for dst parent
      // need a forward reference for it
      sort_ = compute_sort(dst_parent);
      term_dst_parent = ts_.make_inputvar(dst_parent->toString(), sort_);
      w2term_[dst_parent] = term_dst_parent;
    } else {
      term_dst_parent = w2term_.at(dst_parent);
    }

    // expecting bit-vectors, cannot select from a Bool
    assert(term_dst_parent->get_sort()->get_sort_kind() == BV);
    assert(t_->get_sort()->get_sort_kind() == BV);

    tmpterm =
        solver_->make_term(Op(Extract, dst_idx, dst_idx), term_dst_parent);
    Term src_sel_term = solver_->make_term(Op(Extract, src_idx, src_idx), t_);

    // constrain to be equivalent
    ts_.add_constraint(solver_->make_term(Equal, src_sel_term, tmpterm));

  } else {
    tmpterm = t_;
  }

  if (w2term_.find(dst) != w2term_.end()) {
    throw PonoException("CoreIREncoder error. Multiple drivers for "
                        + dst->toString());
  }

  // name and save the value for the dst
  w2term_[dst] = tmpterm;
  ts_.name_term(dst->toString(), tmpterm);
}

Sort CoreIREncoder::compute_sort(CoreIR::Wireable * w)
{
  Type * t = w->getType();
  Sort s;
  if (t->getKind() == CoreIR::Type::TypeKind::TK_Array) {
    // bit-vector sort -- array of bits
    s = solver_->make_sort(BV, t->getSize());
  } else {
    // boolean sort
    s = boolsort_;
  }
  return s;
}

}  // namespace pono
