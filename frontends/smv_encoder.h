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

namespace pono {

class SMVEncoder
{
 public:
  SMVEncoder(std::string filename, pono::RelationalTransitionSystem & rts)
      : rts_(rts), solver_(rts.solver())
  {
    parse(filename);
    processCase();
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
  pono::RelationalTransitionSystem & rts_;
  std::unordered_map<std::string, smt::Term> terms_;
  std::vector<smt::Sort> sortvec_;
  std::vector<smt::Term> propvec_;
  std::unordered_map<std::string, smt::Term> signedbv_;
  std::unordered_map<std::string, smt::Term> unsignedbv_;
  std::deque<std::pair<int, smt::Term>> transterm_;
  std::vector<int> caselist_;
  std::vector<smt::Term> casecheck_;
  std::vector<smt::Term> casestore_;
  // std::unordered_map<int, smt::Term> casecheck_;
  // std::unordered_map<int, smt::Term> casestore_;
  std::vector<std::pair<smt::Term, smt::Term>> caseterm_;
};  // class SMVEncoder
}  // namespace pono
