/*********************                                                  */
/*! \file implicit_predicate_abstractor.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Implicit predicate abstraction based on
**        Abstract Model Checking without Computing the Abstraction
**        Stefano Tonetta
**
**/

#include "modifiers/implicit_predicate_abstractor.h"

#include "assert.h"
#include "utils/logger.h"
#include "smt/available_solvers.h"

using namespace smt;
using namespace std;

namespace pono {

ImplicitPredicateAbstractor::ImplicitPredicateAbstractor(
    const TransitionSystem & conc_ts, TransitionSystem & abs_ts,
    Unroller & un)
    : Abstractor(conc_ts, abs_ts),
      solver_(abs_ts.solver()),
      unroller_(un),
      reducer_(create_solver(solver_->get_solver_enum())),
      abs_rts_(static_cast<RelationalTransitionSystem &>(abs_ts_))
{
  if (conc_ts_.solver() != abs_ts_.solver()) {
    throw PonoException(
        "For simplicity, expecting concrete and abstract system to use same "
        "solver.");
  }

  // TODO: fix abstraction interface
  //       kind of strange to have to pass an empty abstract system

  if (abs_ts_.is_functional()) {
    throw PonoException(
        "Implicit predicate abstraction needs a relational abstract system");
  }

  // assume abstract transition starts empty
  // need to add all state variables and set behavior
  for (auto v : conc_ts_.statevars()) {
    abs_rts_.add_statevar(v, conc_ts_.next(v));
  }
  for (auto v : conc_ts_.inputvars()) {
    abs_rts_.add_inputvar(v);
  }
  // should start with the exact same behavior
  abs_rts_.set_behavior(conc_ts_.init(), conc_ts_.trans());

  do_abstraction();

}

Term ImplicitPredicateAbstractor::abstract(Term & t)
{
  return solver_->substitute(t, abstraction_cache_);
}

Term ImplicitPredicateAbstractor::concrete(Term & t)
{
  return solver_->substitute(t, concretization_cache_);
}

// TODO: somewhere should add predicates from init / prop by default
Term ImplicitPredicateAbstractor::add_predicate(const Term & pred)
{
  assert(abs_ts_.only_curr(pred));
  predicates_.push_back(pred);

  Term rel = predicate_refinement(pred);
  abs_rts_.constrain_trans(rel);
  return rel;
}

bool ImplicitPredicateAbstractor::reduce_predicates(const TermVec & cex,
                                                    const TermVec & new_preds,
                                                    TermVec & out)
{
  assert(new_preds.size());
  Term formula = solver_->make_term(true);

  for (size_t i = 0; i < cex.size(); ++i) {
    formula = solver_->make_term(And, formula, unroller_.at_time(cex[i], i));
    if (i != cex.size() - 1) {
      formula = solver_->make_term(And, formula,
                                   unroller_.at_time(abs_ts_.trans(), i));
    }
  }

  TermVec assumps, red_assumps;
  UnorderedTermMap assumps_to_pred;
  for (const auto &p : new_preds) {
    Term pred_ref = predicate_refinement(p);

    Term a;
    for (size_t i = 0; i+1 < cex.size(); ++i) {
      Term p_i = unroller_.at_time(pred_ref, i);
      a = i==0 ? p_i : solver_->make_term(And, a, p_i);
    }

    assumps.push_back(a);
    assumps_to_pred[a] = p;
  }

  size_t n = out.size();
  reducer_.reduce_assump_unsatcore(formula, assumps, red_assumps, nullptr, 1);
  for (const auto &a : red_assumps) {
    out.push_back(assumps_to_pred.at(a));
  }

  return out.size() > n;
}


void ImplicitPredicateAbstractor::do_abstraction()
{
  logger.log(1, "Generating implicit predicate abstraction.");

  // assume abs_ts_ is relational -- required for this abstraction
  // Note: abs_rts_ is abs_ts_ with a static cast to RelationalTransitionSystem&
  assert(!abs_ts_.is_functional());

  // assume abs_rts_ is a perfect copy currently
  assert(abs_rts_.init() == conc_ts_.init());
  assert(abs_rts_.trans() == conc_ts_.trans());
  assert(abs_rts_.statevars().size() == conc_ts_.statevars().size());
  assert(abs_rts_.inputvars().size() == conc_ts_.inputvars().size());

  Sort boolsort_ = solver_->make_sort(BOOL);

  // create abstract variables for each next state variable
  for (auto sv : conc_ts_.statevars()) {
    if (sv->get_sort() == boolsort_)
    {
      // don't abstract boolean variables
      // but they're implicitly considered predicates
      // which are precise instead of being abstracted
      // so there doesn't need to be a relation added, e.g.
      // P(X') <-> P(X^) is not needed for boolean variables
      assert(abs_ts_.is_curr_var(sv));
      predicates_.push_back(sv);
      continue;
    }
    Term nv = conc_ts_.next(sv);
    // note: this is not a state variable -- using input variable so there's no
    // next
    Term abs_nv = abs_rts_.make_inputvar(nv->to_string() + "^", nv->get_sort());
    // map next var to this abstracted next var
    update_term_cache(nv, abs_nv);
  }

  // TODO: fix the population.
  // Right now state_updates, constraints, and named_terms are not updated
  Term trans = conc_ts_.trans();
  abs_rts_.set_trans(abstract(trans));
  logger.log(5, "Set abstract transition relation to {}", abs_rts_.trans());
}

Term ImplicitPredicateAbstractor::predicate_refinement(const Term & pred)
{
  Term next_pred = abs_ts_.next(pred);
  // constrain next state vars and abstract vars to agree on this predicate
  return solver_->make_term(Equal, next_pred, abstract(next_pred));
}

}  // namespace pono
