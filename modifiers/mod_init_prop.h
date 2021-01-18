/*********************                                                  */
/*! \file mod_init_prop.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Replace initial state and property constraints with boolean variables.
**
**
**/

#pragma once

#include "core/ts.h"
#include "smt-switch/utils.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"

namespace pono {

smt::Term modify_init_and_prop(TransitionSystem & ts, const smt::Term & prop)
{
  logger.log(1, "Modifying init and prop");

  // copy constraints from before we start modifying the system
  std::vector<std::pair<smt::Term, bool>> constraints = ts.constraints();

  // replace prop if it's not already a literal
  smt::Sort boolsort = ts.make_sort(smt::BOOL);
  // NOTE: have to make a new property even if the previous one
  //       was a symbol as well
  //       since we're modifying init, we can't check the property
  //       in init, we have to take a transition first
  smt::Term new_prop = ts.make_statevar("__propvar", boolsort);
  ts.assign_next(new_prop, prop);

  // replace initial states with a single symbol
  smt::Term initstate =
      ts.make_statevar("__initstate", ts.make_sort(smt::BOOL));

  smt::Term init = ts.init();
  smt::TermVec init_constraints;
  conjunctive_partition(init, init_constraints, true);

  // add initial state constraints for initstate
  for (const auto & ic : init_constraints) {
    ts.add_constraint(ts.make_term(smt::Implies, initstate, ic), true);
  }

  for (const auto & e : constraints) {
    ts.add_constraint(ts.make_term(smt::Implies, initstate, e.first), e.second);
  }

  ts.assign_next(initstate, ts.make_term(false));

  // adding the constraints above might have put constraints in init
  // overwrite that now
  ts.set_init(initstate);
  // new_prop is delayed, need to assume it in the initial state
  ts.constrain_init(new_prop);

  return new_prop;
}

// optimization to assume the property in the pre-state
// although confusing, this is sound as long as you always check
// for a property violation over the next-state variables
void prop_in_trans(TransitionSystem & ts, const smt::Term & prop)
{
  // NOTE: CRUCIAL that we pass false here
  // cannot add to init or the next states
  // passing false prevents that
  ts.add_constraint(prop, false);
}

}  // namespace pono
