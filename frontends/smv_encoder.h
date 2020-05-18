#pragma once

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <deque>
#include <future>         // std::async, std::future
#include <chrono>         // std::chrono::milliseconds
#include <unordered_map>
#include "assert.h"
#include "exceptions.h"
#include "rts.h"
#include "smt-switch/smt.h"
#include "smvparser.h"
#include "smvscanner.h"

// #include "bmc.h"
// #include "bmc_simplepath.h"
// #include "defaults.h"
// #include "interpolant.h"
// #include "kinduction.h"
// #include "prop.h"

namespace cosa{
class SMVEncoder
{
 public:
  SMVEncoder(std::string filename, cosa::RelationalTransitionSystem & rts)
      : rts_(rts), solver_(rts.solver())
  {
    parse(filename);
    preprocess();
  };

 public:
  // Important members
  int parse(std::string filename);
  int parseString(std::string newline);
  location loc; 
  void processCase();
  void preprocess();
  smt::TermVec propvec() { return propvec_; }

  smt::SmtSolver & solver_;
  cosa::RelationalTransitionSystem & rts_;
  std::unordered_map<std::string, smt::Sort> sorts_;
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
  void build_case_node(smt::Term cond, smt::Term expr,int lineno);
};  // class SMVEncoder
}  // namespace cosa
