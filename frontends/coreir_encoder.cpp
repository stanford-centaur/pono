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

      // create a state variable for this register
      sort_ = solver_->make_sort(
          BV,
          ipair.second->getModuleRef()->getGenArgs().at("width")->get<int>());
      t_ = ts_.make_state(ipair.first, sort_);
      w2term_[ipair.second] = t_;
      w2term_[ipair.second->sel("out")] = t_;
      ts_.name_term(ipair.second->sel("out")->toString(), t_);
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
              && !instance_of(parent_inst, "coreir", "reg")) {
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
  // TODO: Implement this!
  cout << "got instance " << inst->toString()
       << " but don't know what to do yet!" << endl;
}

}  // namespace cosa
