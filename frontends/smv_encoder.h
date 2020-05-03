#pragma once

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <string>
#include <unordered_map>
#include "assert.h"
#include "exceptions.h"
#include "rts.h"
#include "smt-switch/smt.h"
#include "smvparser.h"
#include "smvscanner.h"

namespace cosa{
class SMVEncoder
{
 public:
  SMVEncoder(std::string filename, cosa::RelationalTransitionSystem & rts)
      : rts_(rts), solver_(rts.solver())
  {
    parse(filename);
  };

 public:
  // Important members
  int parse(std::string filename);
  int parseString(std::string newline);
  void processCase();
  smt::TermVec propvec() { return propvec_; }

  smt::SmtSolver & solver_;
  cosa::RelationalTransitionSystem & rts_;
  std::unordered_map<std::string, smt::Sort> sorts_;
  std::unordered_map<std::string, smt::Term> terms_;
  std::vector<smt::Sort> sortvec_;
  std::vector<smt::Term> propvec_;
  std::unordered_map<std::string, smt::Term>  signedbv_;
  std::unordered_map<std::string, smt::Term>  unsignedbv_;
  std::unordered_map<smt::Term, smt::Term> caseterm_;
};  // class SMVEncoder
}  // namespace cosa
