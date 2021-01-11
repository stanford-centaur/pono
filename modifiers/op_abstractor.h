/*********************                                                  */
/*! \file op_abstractor.h
** \verbatim
** Top contributors (to current version):
**   Hongce Zhang
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract complex operators with uninterpreted functions
**
**
**
**/

#pragma once

#include "smt-switch/utils.h"
#include "engines/ic3base.h"
#include "abstractor.h"

namespace pono {

struct OpAbstract{
  smt::Op op;
  smt::Term result;
  smt::TermVec args;
  
  smt::Term original;
  // statistics
  unsigned refine_count;

  OpAbstract() : refine_count(0) {}
};

class OpAbstractor : public Abstractor
{
public:
  typedef std::unordered_set<smt::PrimOp> OpSet;

  OpAbstractor(const TransitionSystem & conc_ts,
    TransitionSystem & abs_ts,
    const OpSet & op_to_abstract,
    const smt::Term & prop, // it is okay to use bad
    int verbosity);
  
  // return true if it can refine
  // otherwise return false
  bool refine_with_constraints(
    // TODO: how to put the trace in here?
    const ProofGoal * goal_at_init,
    smt::TermVec & out);

protected:
  // find mul/div and replace them with input/output
  void abstract_ts();

  std::vector<OpAbstract> op_abstracted;

}; // class OpAbstractor

} // namespace pono
