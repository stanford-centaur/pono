/*********************                                                  */
/*! \file mod_init_prop.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
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

#include "core/ts.h"
#include "smt-switch/utils.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"

using namespace std;
using namespace smt;

namespace pono {

Term modify_init_and_prop(TransitionSystem & ts, const Term & prop)
{
  logger.log(1, "Modifying init and prop");

  // copy constraints from before we start modifying the system
  TermVec constraints = ts.constraints();

  // replace prop if it's not already a literal
  Sort boolsort = ts.make_sort(BOOL);
  Term new_prop = prop;
  if (!is_lit(prop, boolsort)) {
    new_prop = ts.make_statevar("__propvar", boolsort);
    ts.assign_next(new_prop, prop);
  }

  // replace initial states
  Term initstate1 = ts.make_statevar("__initstate1", ts.make_sort(BOOL));

  Term init = ts.init();
  TermVec init_constraints;
  conjunctive_partition(init, init_constraints, true);

  ts.set_init(initstate1);

  // NOTE: relies on feature of ts to not add constraint to init
  for (const auto & c : constraints) {
    // TODO possibly refactor constraints so next state versions aren't
    // automatically added
    if (ts.no_next(c)) {
      ts.add_constraint(ts.make_term(Implies, initstate1, c), false);
    }
  }

  ts.assign_next(initstate1, ts.make_term(false));
  ts.constrain_init(new_prop);

  // add initial state constraints for initstate1
  for (const auto & ic : init_constraints) {
    ts.add_constraint(ts.make_term(Implies, initstate1, ic), false);
  }

  return new_prop;
}

}  // namespace pono
