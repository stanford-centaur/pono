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

#include <cassert>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>

#include "core/ts.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"

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
  const std::vector<smt::TermVec> & justicevec() const { return justicevec_; };
  const smt::TermVec & fairvec() const { return fairvec_; };
  const smt::TermVec & inputsvec() const { return inputsvec_; }
  const smt::TermVec & statesvec() const { return statesvec_; }
  const std::map<uint64_t, smt::Term> & no_next_statevars() const
  {
    return no_next_states_;
  }
  const std::unordered_map<std::string, std::string> & get_symbol_map() const
  {
    return symbol_map_;
  }
  const smt::UnorderedTermSet & initialized_statevars() const
  {
    return initialized_states_;
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
  // record the renaming done by the `preprocess` pass
  std::unordered_map<uint64_t, std::string> state_renaming_table_;
  // a mapping from the internally-assigned names ("intput/state{btor2_id}")
  // to the original names in Btor2
  std::unordered_map<std::string, std::string> symbol_map_;
  smt::UnorderedTermSet initialized_states_;

  std::unordered_map<int, smt::Sort> sorts_;
  std::unordered_map<int, smt::Term> terms_;

  // properties, justice, and fairness constraints
  smt::TermVec propvec_;
  std::vector<smt::TermVec> justicevec_;
  smt::TermVec fairvec_;
};
}  // namespace pono
