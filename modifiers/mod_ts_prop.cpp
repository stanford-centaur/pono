/*********************                                                  */
/*! \file mod_ts_prop.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Functions for modifying ts and property
**
**
**/

#include "modifiers/mod_ts_prop.h"

#include "core/rts.h"
#include "smt-switch/utils.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"
#include "utils/ts_manipulation.h"

using namespace smt;
using namespace std;

namespace pono {

TransitionSystem pseudo_init_and_prop(TransitionSystem & ts, Term & prop)
{
  logger.log(1, "Modifying init and prop");

  RelationalTransitionSystem rts =
      ts.is_functional() ? RelationalTransitionSystem(ts.solver())
                         : static_cast<RelationalTransitionSystem &>(ts);

  // make it relational if it wasn't already
  if (ts.is_functional()) {
    for (const auto & sv : ts.statevars()) {
      rts.add_statevar(sv, ts.next(sv));
    }
    for (const auto & iv : ts.inputvars()) {
      rts.add_inputvar(iv);
    }
    rts.set_init(ts.init());
    rts.set_trans(ts.trans());
  }

  assert(!rts.is_functional());

  Sort boolsort = rts.make_sort(BOOL);

  // Save original initial state constraints
  TermVec init_conjuncts;
  conjunctive_partition(rts.init(), init_conjuncts, true);

  // create a pseudo initial state
  Term pseudo_init = rts.make_statevar("__pseudo_init", boolsort);
  Term not_pseudo_init = rts.make_term(Not, pseudo_init);

  // guard property with it
  prop = rts.make_term(Implies, not_pseudo_init, prop);

  // create a fake transition from pseudo initial state
  // and enforce that the second state is constrained by the original
  // initial state constraints
  // NOTE: need to guard transition relation as well for correctness
  TermVec trans_conjuncts;
  conjunctive_partition(rts.trans(), trans_conjuncts, true);
  // reset transition relation
  rts.set_trans(rts.make_term(true));
  // add fake transition from pseudo_init to actual initial states
  for (const auto & ic : init_conjuncts) {
    assert(rts.only_curr(ic));
    rts.constrain_trans(rts.make_term(Implies, pseudo_init, rts.next(ic)));
  }
  // add regular transitions back in
  for (const auto & tc : trans_conjuncts) {
    rts.constrain_trans(rts.make_term(Implies, not_pseudo_init, tc));
  }

  rts.set_init(pseudo_init);
  rts.assign_next(pseudo_init, rts.make_term(false));

  // now create a property monitor
  Term new_prop = rts.make_statevar("__prop_monitor", boolsort);
  if (rts.only_curr(prop)) {
    rts.add_invar(rts.make_term(Equal, new_prop, prop));
  } else {
    // property monitor
    rts.assign_next(new_prop, prop);
  }
  // need to assume in the pseudo initial state
  // or else there's a trivial counterexample
  rts.constrain_init(new_prop);

  // update the property in-place
  std::swap(prop, new_prop);

  return rts;
}

void prop_in_trans(TransitionSystem & ts, const Term & prop)
{
  // NOTE: CRUCIAL that we pass false here
  // cannot add to init or the next states
  // passing false prevents that
  ts.add_constraint(prop, false);
}

TransitionSystem promote_inputvars(const TransitionSystem & ts)
{
  SmtSolver solver = ts.solver();
  TransitionSystem new_ts = create_fresh_ts(ts.is_functional(), solver);
  assert(new_ts.is_functional() == ts.is_functional());

  // copy over all state variables
  for (const auto & sv : ts.statevars()) {
    new_ts.add_statevar(sv, ts.next(sv));
  }

  // copy over inputs but make them statevars
  for (const auto & iv : ts.inputvars()) {
    new_ts.add_inputvar(iv);
    new_ts.promote_inputvar(iv);
  }

  // set init
  new_ts.set_init(ts.init());

  // copy over state updates
  for (const auto & elem : ts.state_updates()) {
    new_ts.assign_next(elem.first, elem.second);
  }

  // need to re-evaluate all constraints that used to be over inputs
  for (const auto & elem : ts.constraints()) {
    new_ts.add_constraint(elem.first, elem.second);
  }

  // relational systems could have things added by constrain_trans
  if (!new_ts.is_functional()) {
    RelationalTransitionSystem & rts_view =
        static_cast<RelationalTransitionSystem &>(new_ts);
    rts_view.set_trans(ts.trans());
  }

  assert(!new_ts.inputvars().size());
  return new_ts;
}

}  // namespace pono
