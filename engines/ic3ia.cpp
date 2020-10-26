/*********************                                                  */
/*! \file ic3ia.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief IC3 via Implicit Predicate Abstraction (IC3IA) implementation
**        based on
**
**        IC3 Modulo Theories via Implicit Predicate Abstraction
**            -- Alessandro Cimatti, Alberto Griggio,
**               Sergio Mover, Stefano Tonetta
**
**        and the open source implementation:
**
**        https://es-static.fbk.eu/people/griggio/ic3ia/index.html
**/

#include "engines/ic3ia.h"

#include "utils/term_analysis.h"

using namespace smt;
using namespace std;

namespace pono {

IC3IA::IC3IA(Property & p, SolverEnum se)
    : super(p, se), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

IC3IA::IC3IA(Property & p, const SmtSolver & slv)
    : super(p, slv), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

IC3IA::IC3IA(const PonoOptions & opt, Property & p, const SolverEnum se)
    : super(opt, p, se), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

IC3IA::IC3IA(const PonoOptions & opt, Property & p, const SmtSolver & slv)
    : super(opt, p, slv), abs_ts_(ts_.solver()), ia_(ts_, abs_ts_)
{
  initialize();
}

// protected methods
void IC3IA::initialize()
{
  boolsort_ = solver_->make_sort(BOOL);
  // add all the predicates from init and property
  UnorderedTermSet preds;
  get_predicates(ts_.init(), boolsort_, preds, true);
  get_predicates(bad_, boolsort_, preds, true);
  for (auto p : preds) {
    add_predicate(p);
  }

  // TODO: set semantics for trans, etc...
  throw PonoException("NYI");
}

void IC3IA::add_predicate(const Term & pred)
{
  assert(abs_ts_.only_curr(pred));
  // add predicate to abstraction and get the new constraint
  Term predabs_rel = ia_.add_predicate(pred);
  // refine the transition relation incrementally
  // by adding a new constraint
  assert(!solver_context_);  // should be at context 0
  solver_->assert_formula(
      solver_->make_term(Implies, trans_label_, predabs_rel));

  // create a new state variable for the IC3 procedure
  // used to represent the predicate
  // NOTE: don't need to actually add to the TS
  //       kept separate for predicate abstraction
  // TODO: add infrastructure for getting a fresh symbol -- might fail if name
  // already exists
  string name = "__pred" + std::to_string(pred_statevar.size());
  Term pred_state = solver_->make_symbol(name, boolsort_);
  Term pred_next = solver_->make_symbol(name + ".next", boolsort_);
  pred_statevars_.push_back(pred_state);
  predstate2next_[pred_state] = pred_next;

  // associate the fresh state vars with predicate applied to current / next
  // states
  solver_->assert_formula(solver_->make_term(Implies, pred_state, pred));
  // NOTE: this predicate is over the original next state vars (not abstracted)
  solver_->assert_formula(
      solver_->make_term(Implies, pred_state, ts_.next(pred)));
}

}  // namespace pono
