#pragma once

#include <stdio.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include "assert.h"
#include "exceptions.h"
#include "rts.h"
#include "smt-switch/smt.h"
#include "smvparser.h"

class smvEncoder
{
 public:
  smvEncoder(std::string filename, cosa::RelationalTransitionSystem & rts)
      : rts_(rts), solver_(rts.solver())
  {
    parse(filename);
  };

 public:
  // Important members
  // void read(std::string filename);
  void parse(std::string filename);

  smt::TermVec propvec() { return propvec_; }

  smt::SmtSolver & solver_;
  cosa::RelationalTransitionSystem & rts_;
  std::unordered_map<std::string, smt::Sort> sorts_;
  std::unordered_map<std::string, smt::Term> terms_;
  std::vector<smt::Sort> sortvec_;
  std::vector<smt::Term> propvec_;

  int id;
  // std::string symbol_;  // namespace cosa
};
