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
  smt::Term new_prop = prop;
  if (!is_lit(prop, boolsort)) {
    new_prop = ts.make_statevar("__propvar", boolsort);
    ts.assign_next(new_prop, prop);
  }

  // replace initial states
  smt::Term initstate1 = ts.make_statevar("__initstate1",
                                          ts.make_sort(smt::BOOL));

  smt::Term init = ts.init();
  smt::TermVec init_constraints;
  conjunctive_partition(init, init_constraints, true);

  // NOTE: relies on feature of ts to not add constraint to init
  for (const auto & e : constraints) {
    ts.add_constraint(ts.make_term(smt::Implies, initstate1, e.first),
                      e.second);
  }

  ts.assign_next(initstate1, ts.make_term(false));

  // adding the constraints above might have put constraints in init
  // overwrite that now
  ts.set_init(initstate1);
  ts.constrain_init(new_prop);

  // add initial state constraints for initstate1
  for (const auto & ic : init_constraints) {
    ts.add_constraint(ts.make_term(smt::Implies, initstate1, ic), false);
  }

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
