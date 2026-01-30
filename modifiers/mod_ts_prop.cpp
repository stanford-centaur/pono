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

Term pseudo_init_and_prop(TransitionSystem & ts, const Term & prop)
{
  assert(!ts.is_functional());

  logger.log(1, "Modifying init and prop");
  RelationalTransitionSystem & rts =
      static_cast<RelationalTransitionSystem &>(ts);
  Sort boolsort = rts.make_sort(BOOL);

  // Save original initial state constraints
  TermVec init_conjuncts;
  conjunctive_partition(rts.init(), init_conjuncts, true);

  // create a pseudo initial state
  Term pseudo_init = rts.make_statevar("__pseudo_init", boolsort);
  Term not_pseudo_init = rts.make_term(Not, pseudo_init);

  // guard property with it
  Term guarded_prop = rts.make_term(Implies, not_pseudo_init, prop);

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
  if (rts.only_curr(guarded_prop)) {
    rts.add_invar(rts.make_term(Equal, new_prop, guarded_prop));
  } else {
    // property monitor
    rts.assign_next(new_prop, guarded_prop);
  }
  // need to assume in the pseudo initial state
  // or else there's a trivial counterexample
  rts.constrain_init(new_prop);

  return new_prop;
}

void prop_in_trans(TransitionSystem & ts, const Term & prop)
{
  // NOTE: CRUCIAL that we pass false here
  // cannot add to init or the next states
  // passing false prevents that
  ts.add_constraint(prop, false);
}

TransitionSystem promote_inputvars(const TransitionSystem & ts,
                                   const UnorderedTermSet & ivs_to_promote)
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
    if (ivs_to_promote.find(iv) != ivs_to_promote.end()) {
      new_ts.promote_inputvar(iv);
    }
  }

  // set init
  new_ts.set_init(ts.init());

  // copy over state updates
  for (const auto & elem : ts.state_updates()) {
    new_ts.assign_next(elem.first, elem.second);
  }

  // relational systems could have things added by constrain_trans;
  // this step has to be performed before re-adding constraints as new
  // constraints might be added to trans()
  if (!new_ts.is_functional()) {
    RelationalTransitionSystem & rts_view =
        static_cast<RelationalTransitionSystem &>(new_ts);
    rts_view.set_trans(ts.trans());
  }

  // need to re-evaluate all constraints that used to be over inputs
  for (const auto & elem : ts.constraints()) {
    new_ts.add_constraint(elem.first, elem.second);
  }

  return new_ts;
}

TransitionSystem promote_inputvars(const TransitionSystem & ts)
{
  TransitionSystem ret = promote_inputvars(ts, ts.inputvars());
  assert(ret.inputvars().empty());
  return ret;
}

void make_trans_total(TransitionSystem & ts, Term & prop)
{
  if (ts.is_functional() && ts.constraints().empty()) {
    // already right-total
    return;
  }

  // create a new ts without constraints
  TransitionSystem new_ts = create_fresh_ts(ts.is_functional(), ts.solver());

  // copy all state and input variables
  for (const auto & sv : ts.statevars()) {
    new_ts.add_statevar(sv, ts.next(sv));
  }
  for (const auto & iv : ts.inputvars()) {
    new_ts.add_inputvar(iv);
  }
  // set init
  new_ts.set_init(ts.init());

  // instrument ts with a new state var valid
  Term valid;
  std::size_t id = 0;
  while (true) {
    try {
      valid = new_ts.make_statevar("__pono_valid_trans_" + std::to_string(id),
                                   ts.make_sort(smt::BOOL));
      break;
    }
    catch (SmtException & e) {
      ++id;
    }
  }
  new_ts.constrain_init(valid);  // valid@0 = true

  if (ts.is_functional()) {
    for (const auto & elem : ts.state_updates()) {
      new_ts.assign_next(elem.first, elem.second);
    }

    vector<Term> constraints;
    constraints.reserve(ts.constraints().size());
    for (const auto & elem : ts.constraints()) {
      constraints.push_back(elem.first);
    }
    Term valid_next = new_ts.make_term(
        And,
        valid,
        constraints.size() == 1 ? constraints.at(0)
                                : new_ts.make_term(And, constraints));
    new_ts.assign_next(valid, valid_next);
    // `valid` is later used for modifying `prop`
    // the constraints at current time step have to be included
    valid = valid_next;
  } else {
    Term new_trans = new_ts.make_term(
        Or,
        new_ts.make_term(
            And,
            valid,
            new_ts.make_term(Equal, new_ts.next(valid), ts.trans())),
        new_ts.make_term(And,
                         new_ts.make_term(Not, valid),
                         new_ts.make_term(Not, new_ts.next(valid))));
    static_cast<RelationalTransitionSystem &>(new_ts).set_trans(new_trans);
  }
  assert(new_ts.is_right_total());

  // update ts and prop
  ts = new_ts;
  prop = ts.make_term(Implies, valid, prop);
}

}  // namespace pono
