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

namespace cosa {

// static functions

Module * CoreIREncoder::read_coreir_file(Context * c, std::string filename)
{
  Module * m;
  // load file
  if (!loadFromFile(c, filename, &m)) {
    c->die();
    throw CosaException("Error reading CoreIR file: " + filename);
  }
  return m;
}

// member functions

void CoreIREncoder::parse(std::string filename)
{
  top_ = read_coreir_file(c_, filename);

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

  // create inputs for interface inputs and states for clocks
  for (auto elem : def_->getInterface()->getSelects()) {
    if (elem.second->getType()->toString() == "coreir.clk") {
      t_ = ts_.make_state(elem.first, solver_->make_sort(BOOL));
      w2term_[elem.second] = t_;
      num_clocks_++;
    }
    type_ = elem.second->getType();
    if (type_->isInput() || type_->isInOut()) {
      t_ = ts_.make_input(elem.first, solver_->make_sort(BV, type_->getSize()));
      w2term_[elem.second] = t_;
    }
  }

  // can't abstract the clock if there's more than one
  can_abstract_clock_ &= (num_clocks_ <= 1);

  // create registers and store the number of inputs for each instance
  vector<Instance *> instances;
  unordered_map<Instance *, size_t> num_inputs;
  bool async;
  for (auto ipair : def_->getInstances()) {
    type_ = ipair.second->getType();
    if (instance_of(ipair.second, "coreir", "reg")
        || (async = instance_of(ipair.second, "coreir", "arst"))) {
      // cannot abstract clock if there are asynchronous resets
      can_abstract_clock_ &= !async;

      // put registers into instances first (processed last)
      instances.push_back(ipair.second);
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

  logger.log(1,
             "INFO {} abstract clock for CoreIR file {}",
             can_abstract_clock_ ? "can" : "cannot",
             filename);

  // process the rest in topological order
  size_t processed_instances = 0;
  unordered_map<Instance *, set<Wireable *>> covered_inputs;
  while (instances.size()) {
    inst_ = instances.back();
    instances.pop_back();

    process_instance(inst_);
    processed_instances++;

    // check everything connected to outputs
    // append to list if all other inputs have been seen already

    Wireable * dst;
    Wireable * parent;
    Instance * parent_inst;
    for (auto conn : inst_->getLocalConnections()) {
      dst = conn.second;
      type_ = dst->getType();
      if (type_->isInput() || type_->isInOut()) {
        // parent is either an instance or a top-level input
        parent = dst->getTopParent();
        if (Instance::classof(parent)) {
          parent_inst = dyn_cast<Instance>(parent);
          // states have already been added so ignore those
          if (!instance_of(parent_inst, "coreir", "reg")
              && !instance_of(parent_inst, "coreir", "reg_arst")) {
            covered_inputs[parent_inst].insert(dst);
            if (num_inputs[parent_inst]
                == covered_inputs.at(parent_inst).size()) {
              instances.push_back(parent_inst);
            }
          }
        }
        // TODO: consider have an else case here to create the global inputs
        // instead of doing it ahead of time earlier
      }
    }
  }

  if (processed_instances != def_->getInstances().size()) {
    throw CosaException("Issue: not all instances processed in CoreIR Encoder");
  }
}

void CoreIREncoder::process_instance(CoreIR::Instance * inst)
{
  mod_ = inst->getModuleRef();
  string nsname = mod_->getGenerator()->getNamespace()->getName();
  string name = mod_->getGenerator()->getName();
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
    // NOTE: inputs to registers are not wired up until later
    sort_ = solver_->make_sort(
        BV, inst->getModuleRef()->getGenArgs().at("width")->get<int>());
    t_ = ts_.make_state(inst->toString(), sort_);
  } else if (nsname == "coreir" && name == "const") {
    size_t w = mod_->getGenArgs().at("width")->get<int>();
    sort_ = solver_->make_sort(BV, w);
    t_ = solver_->make_term(
        (inst->getModArgs().at("value"))->get<BitVec>().binary_string(),
        sort_,
        2);
  } else if (nsname == "corebit" && name == "const") {
    sort_ = solver_->make_sort(BOOL);
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
    t_ = ts_.make_input(inst->toString(), sort_);
  } else if (nsname == "corebit" && name == "undriven") {
    t_ = ts_.make_input(inst->toString(), boolsort_);
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
    return;
  }

  if (!t_) {
    throw CosaException("CoreIREncoder error: no term created for module type: "
                        + mod_->getName());
  }

  w2term_[inst] = t_;
  w2term_[inst->sel("out")] = t_;
  ts_.name_term(inst->sel("out")->toString(), t_);

  Term tmpterm;
  for (auto conn : inst->sel("out")->getLocalConnections()) {
    if (isa<CoreIR::Select>(conn.first)
        && isNumber((sel_ = cast<CoreIR::Select>(conn.first))->getSelStr())) {
      size_t idx = stoi(sel_->getSelStr());
      tmpterm = solver_->make_term(Op(Extract, idx, idx), t_);
      tmpterm = solver_->make_term(Equal, tmpterm, bv1_);
      w2term_[conn.first] = tmpterm;
    } else {
      tmpterm = t_;
    }

    // name and save the value for the input
    w2term_[conn.second] = tmpterm;
    ts_.name_term(conn.second->toString(), tmpterm);
  }
}

}  // namespace cosa
