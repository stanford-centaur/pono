/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#pragma once


extern "C" {
#include "btor2parser/btor2parser.h"
}

#include <stdio.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include "assert.h"

#include "core/ts.h"
#include "utils/exceptions.h"

#include "smt-switch/smt.h"

namespace pono {
class BTOR2Encoder
{
 public:
  BTOR2Encoder(std::string filename, TransitionSystem & ts)
      : ts_(ts), solver_(ts.solver())
  {
    preprocess(filename);
    parse(filename);
  };

  const smt::TermVec & propvec() const { return propvec_; };
  const smt::TermVec & justicevec() const { return justicevec_; };
  const smt::TermVec & fairvec() const { return fairvec_; };
  const smt::TermVec & inputsvec() const { return inputsvec_; }
  const smt::TermVec & statesvec() const { return statesvec_; }
  const std::map<uint64_t, smt::Term> & no_next_statevars() const
  {
    return no_next_states_;
  }

 protected:
  // converts booleans to bitvector of size one
  smt::Term bool_to_bv(const smt::Term & t) const;
  // converts bitvector of size one to boolean
  smt::Term bv_to_bool(const smt::Term & t) const;
  // lazy conversion
  // takes a list of booleans / bitvectors of size one
  // and lazily converts them to the majority
  smt::TermVec lazy_convert(const smt::TermVec &) const;
  
  // preprocess a btor2 file
  void preprocess(const std::string & filename);
  // parse a btor2 file
  void parse(const std::string filename);

  // Important members
  const smt::SmtSolver & solver_;
  pono::TransitionSystem & ts_;

  // vectors of inputs and states
  // maintains the order from the btor file
  smt::TermVec inputsvec_;
  smt::TermVec statesvec_;
  std::map<uint64_t, smt::Term> no_next_states_;
  std::unordered_map<uint64_t, std::string> state_renaming_table;

  // Useful variables
  smt::Sort linesort_;
  smt::TermVec termargs_;
  std::unordered_map<int, smt::Sort> sorts_;
  std::unordered_map<int, smt::Term> terms_;
  std::string symbol_;

  smt::TermVec propvec_;
  smt::TermVec justicevec_;
  smt::TermVec fairvec_;

  Btor2Parser * reader_;
  Btor2LineIterator it_;
  Btor2Line * l_;
  size_t i_;
  int64_t idx_;
  bool negated_;
  size_t witness_id_{ 0 };  ///< id of any introduced witnesses for properties
};
}  // namespace pono
