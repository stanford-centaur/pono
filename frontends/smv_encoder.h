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
#include "exceptions.h"
#include "rts.h"
#include "smt-switch/smt.h"
#include "smvparser.h"
#include "smvscanner.h"

namespace pono{
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
    preprocess();
    processCase();
  };

 public:
  // Important members
  int parse(std::string filename);
  int parse_flat(std::istream& s);
  int parseString(std::string newline);
  location loc;
  void processCase();
  void preprocess();
  smt::TermVec propvec() { return propvec_; }

  smt::SmtSolver & solver_;
  pono::RelationalTransitionSystem & rts_;
  std::unordered_map<std::string, smt::Term> terms_;
  std::vector<smt::Sort> sortvec_;
  std::vector<smt::Term> propvec_;
  std::unordered_map<std::string, smt::Term>  signedbv_;
  std::unordered_map<std::string, smt::Term>  unsignedbv_;
  std::deque<std::pair<int, smt::Term>> transterm_;
  std::vector<int> caselist_;
  std::unordered_map<int, smt::Term> casecheck_;
  std::unordered_map<int, smt::Term> casestore_;
  std::vector<std::pair<smt::Term, smt::Term>> caseterm_;
  std::unordered_map<std::string,module_node*> module_list;
  bool module_flat;

  std::vector<cosa::SMVnode*> define_list_;
  std::vector<cosa::SMVnode*> assign_list_;
  std::vector<cosa::SMVnode*> ivar_list_;
  std::vector<cosa::var_node_c*> var_list_;
  std::vector<cosa::SMVnode*> frozenvar_list_;
  std::vector<cosa::SMVnode*> init_list_;
  std::vector<cosa::SMVnode*> trans_list_;
  std::vector<cosa::SMVnode*> invar_list_;
  std::vector<cosa::SMVnode*> invarspec_list_;
};  // class SMVEncoder
}  // namespace pono
