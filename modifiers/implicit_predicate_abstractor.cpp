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

bool contains_var(const Term & p, const UnorderedTermSet & vars)
{
  UnorderedTermSet free_vars;
  get_free_symbolic_consts(p, free_vars);
  for (const auto & v : free_vars)
  {
    if (vars.find(v) != vars.end())
    {
      return true;
    }
  }
  return false;
}

ImplicitPredicateAbstractor::ImplicitPredicateAbstractor(
    const TransitionSystem & conc_ts, TransitionSystem & abs_ts, Unroller & un)
    : Abstractor(conc_ts, abs_ts),
      solver_(abs_ts.solver()),
      unroller_(un),
      reducer_(create_solver(solver_->get_solver_enum(), false, true, false)),
      to_reducer_(reducer_),
      abs_rts_(static_cast<RelationalTransitionSystem &>(abs_ts_)),
      red_can_reset_(true)  // start by assuming it can reset
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

}

Term ImplicitPredicateAbstractor::abstract(Term & t)
{
  assert(abstracted_);
  return solver_->substitute(t, abstraction_cache_);
}

Term ImplicitPredicateAbstractor::concrete(Term & t)
{
  assert(abstracted_);
  return solver_->substitute(t, concretization_cache_);
}

Term ImplicitPredicateAbstractor::predicate_refinement(const Term & pred)
{
  assert(abstracted_);
  Term next_pred = abs_ts_.next(pred);
  // constrain next state vars and abstract vars to agree on this predicate
  return solver_->make_term(Equal, next_pred, abstract(next_pred));
}

bool ImplicitPredicateAbstractor::reduce_predicates(const TermVec & cex,
                                                    const TermVec & new_preds,
                                                    TermVec & out)
{
  assert(abstracted_);
  assert(new_preds.size());

  Term formula = solver_->make_term(true);

  for (size_t i = 0; i < cex.size(); ++i) {
    formula = solver_->make_term(And, formula, unroller_.at_time(cex[i], i));
    if (i != cex.size() - 1) {
      formula = solver_->make_term(And, formula,
                                   unroller_.at_time(abs_ts_.trans(), i));
    }
  }

  TermVec assumps;
  for (const auto &p : new_preds) {
    Term pred_ref = predicate_refinement(p);

    Term a;
    for (size_t i = 0; i+1 < cex.size(); ++i) {
      Term p_i = unroller_.at_time(pred_ref, i);
      a = i==0 ? p_i : solver_->make_term(And, a, p_i);
    }

    assumps.push_back(to_reducer_.transfer_term(a, BOOL));
  }
  assert(assumps.size() == new_preds.size());

  size_t n = out.size();
  bool reset = reset_reducer();
  if (!reset) {
    reducer_->push();
  }
  reducer_->assert_formula(to_reducer_.transfer_term(formula));
  Result res = reducer_->check_sat_assuming(assumps);
  if (!res.is_unsat()) {
    // TODO: investigate this in more detail
    // it doesn't seem like this should happen
    out = new_preds;
    return false;
  }

  UnorderedTermSet core;
  try {
    reducer_->get_unsat_assumptions(core);
    for (size_t i = 0; i < assumps.size(); ++i) {
      if (core.find(assumps[i]) != core.end()) {
        out.push_back(new_preds[i]);
      }
      else if (contains_var(new_preds[i], important_vars_))
      {
        logger.log(1, "IA Keeping ImpVar Pred: {}", new_preds[i]);
        out.push_back(new_preds[i]);
      }
    }
  }
  catch (SmtException & e) {
    logger.log(2, "Failed to reduce predicates");
  }

  if (!reset) {
    reducer_->pop();
  }

  return out.size() > n;
}

UnorderedTermSet ImplicitPredicateAbstractor::do_abstraction()
{
  logger.log(1, "Generating implicit predicate abstraction.");

  abstracted_ = true;

  UnorderedTermSet conc_predicates;

  // need to add all state variables and set behavior
  // due to incrementality, most of the variables will be already present
  // so make sure to check that
  const UnorderedTermSet & abs_statevars = abs_rts_.statevars();
  for (const auto &v : conc_ts_.statevars()) {
    if (abs_statevars.find(v) == abs_statevars.end()) {
      abs_rts_.add_statevar(v, conc_ts_.next(v));
    }
  }
  const UnorderedTermSet & abs_inputs = abs_rts_.inputvars();
  for (const auto &v : conc_ts_.inputvars()) {
    if (abs_inputs.find(v) == abs_inputs.end()) {
      abs_rts_.add_inputvar(v);
    }
  }
  // should start with the exact same behavior
  abs_rts_.set_behavior(conc_ts_.init(), conc_ts_.trans());

  // assume abs_ts_ is relational -- required for this abstraction
  // Note: abs_rts_ is abs_ts_ with a static cast to RelationalTransitionSystem&
  assert(!abs_ts_.is_functional());

  // assume abs_rts_ is a perfect copy currently
  assert(abs_rts_.init() == conc_ts_.init());
  assert(abs_rts_.trans() == conc_ts_.trans());
  assert(abs_rts_.statevars().size() == conc_ts_.statevars().size());
  // due to incrementality, there are more input variables in abs_rts
  assert(abs_rts_.inputvars().size() >= conc_ts_.inputvars().size());

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
      conc_predicates.insert(sv);
      continue;
    }

    Term nv = conc_ts_.next(sv);
    // note: this is not a state variable -- using input variable so there's no
    // next

    // for incrementality
    // only create new variables if not present in cache
    if (abstraction_cache_.find(nv) == abstraction_cache_.end()) {
      Term abs_nv = abs_rts_.make_inputvar(nv->to_string() + "^",
                                           nv->get_sort());
      // map next var to this abstracted next var
      update_term_cache(nv, abs_nv);
    }
  }

  // TODO: fix the population.
  // Right now state_updates, constraints, and named_terms are not updated
  Term trans = conc_ts_.trans();
  abs_rts_.set_trans(abstract(trans));
  logger.log(3, "Set abstract transition relation to {}", abs_rts_.trans());

  return conc_predicates;
}

}  // namespace pono
