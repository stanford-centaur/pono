#pragma once

#include <stdio.h>

#include <chrono>  // std::chrono::milliseconds
#include <deque>
#include <fstream>
#include <future>  // std::async, std::future
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "assert.h"
#include "core/rts.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"

#include "frontends/smvscanner.h"
#include "smvparser.h"
#include "smvscanner.h"

namespace pono {

class SMVEncoder
{
 public:
  SMVEncoder(std::string filename, pono::RelationalTransitionSystem & rts)
      : rts_(rts), solver_(rts.solver())
  {
    module_flat = false;
    parse(filename);
    preprocess();
    module_flat = true;
    std::string flatten =filename.substr(0,filename.find_last_of(".")) + "_flatten.txt";
    std::ofstream ofile(flatten);
    ofile << preprocess().str();
    ofile.close();
    processCase();
  };

 public:
  // Important members
  int parse(std::string filename);
  int parse_flat(std::istream& s);
  smt::Term parseString(std::string newline);
  location loc;
  void processCase();
  std::stringstream preprocess();
  smt::TermVec propvec() { return propvec_; }

  smt::Term parse_term;
  const smt::SmtSolver & solver_;
  pono::RelationalTransitionSystem & rts_;
  std::unordered_map<std::string, smt::Term> terms_;
  std::vector<smt::Term> propvec_;
  ///< signedbv_: to store signed bitvector for type checking
  ///< unsignedbv_: to store unsigned bitvector for type checking
  ///< arrayty_: used to store word array nad its corresponding name for 
  std::unordered_map<std::string, smt::Term>  signedbv_;
  std::unordered_map<std::string, smt::Term>  unsignedbv_;
  std::deque<std::pair<int, smt::Term>> transterm_;
  std::unordered_map<std::string, SMVnode::Type>  arrayty_;
  ///< casecheck_: vector of booleans, each element is an Or of all the conditions in a case statement.
  ///< caseterm_: used to temporaily store each statement in case body before future process check that the conditions cover all possibilities (required by nuXmv manual)
  ///< casestore_: used to temporarily store constraint with ite tree, if the condition is satisfied, add the constraint
  std::vector<smt::Term> casecheck_;
  std::vector<smt::Term> casestore_;
  std::vector<std::pair<SMVnode*,SMVnode*>> caseterm_;
  ///< module_list: map from module name to module node
  std::unordered_map<std::string,module_node*> module_list;
  // indicate whether needs to flatten module first
  bool module_flat;

  std::vector<pono::SMVnode*> define_list_;
  std::vector<pono::SMVnode*> assign_list_;
  std::vector<pono::SMVnode*> ivar_list_;
  std::vector<pono::var_node_c*> var_list_;
  std::vector<pono::SMVnode*> frozenvar_list_;
  std::vector<pono::SMVnode*> init_list_;
  std::vector<pono::SMVnode*> trans_list_;
  std::vector<pono::SMVnode*> invar_list_;
  std::vector<pono::SMVnode*> invarspec_list_;
};  // class SMVEncoder
}  // namespace pono
