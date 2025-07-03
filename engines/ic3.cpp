/*********************                                                  */
/*! \file ic3.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Bit-level IC3 implementation using the IC3Base abstract base class
**/

#include "engines/ic3.h"

#include <cassert>

using namespace smt;
using namespace std;

namespace pono {

/** IC3 Implementation */

IC3::IC3(const SafetyProperty & p,
         const TransitionSystem & ts,
         const smt::SmtSolver & s,
         PonoOptions opt)
    : super(p, ts, s, opt)
{
  engine_ = Engine::IC3_BOOL;
}

IC3Formula IC3::get_model_ic3formula() const
{
  // expecting all solving in IC3 to be done at context level > 0
  // so if we're getting a model we should not be at context 0
  assert(solver_context_);

  const UnorderedTermSet & statevars = ts_.statevars();
  TermVec children;
  children.reserve(statevars.size());
  for (const auto & sv : ts_.statevars()) {
    if (solver_->get_value(sv) == solver_true_) {
      children.push_back(sv);
    } else {
      children.push_back(solver_->make_term(Not, sv));
    }
  }

  return ic3formula_conjunction(children);
}

bool IC3::ic3formula_check_valid(const IC3Formula & u) const
{
  // check that children are literals
  Op op;
  for (const auto & c : u.children) {
    if (!is_lit(c, boolsort_)) {
      return false;
    }
  }

  // for now not checking the actual term (e.g. u.term)
  // it's somewhat hard if the underlying solver does rewriting

  // got through all checks without failing
  return true;
}

void IC3::predecessor_generalization(size_t i,
                                     const Term & c,
                                     IC3Formula & pred)
{
  // TODO: change this so we don't have to depend on the solver context to be
  // sat
  assert(i > 0);
  assert(!pred.disjunction);

  const UnorderedTermSet & statevars = ts_.statevars();
  TermVec input_lits = get_input_values();
  TermVec next_lits = get_next_state_values();
  const TermVec & cube_lits = pred.children;

  if (i == 1) {
    // don't need to generalize if i == 1
    // the predecessor is an initial state
    return;
  }

  Term formula = make_and(input_lits);
  if (ts_.is_deterministic()) {
    // NOTE: need to use full trans, not just trans_label_ here
    //       because we are passing it to the reducer_
    formula = solver_->make_term(And, formula, ts_.trans());
    formula =
        solver_->make_term(And, formula, solver_->make_term(Not, ts_.next(c)));
  } else {
    formula = solver_->make_term(And, formula, make_and(next_lits));

    // preimage formula
    // NOTE: because this will be negated, it is important to use the
    // whole frame and trans, not just the labels
    // because label is: trans_label_ -> trans
    // so if it is negated, that doesn't force trans to be false
    // the implication could be more efficient than iff so we want to leave it
    // that way
    Term pre_formula = get_frame_term(i - 1);
    pre_formula = solver_->make_term(And, pre_formula, ts_.trans());
    pre_formula =
        solver_->make_term(And, pre_formula, solver_->make_term(Not, c));
    pre_formula = solver_->make_term(And, pre_formula, ts_.next(c));
    formula =
        solver_->make_term(And, formula, solver_->make_term(Not, pre_formula));
  }
  // TODO: consider adding functional pre-image here
  //       not sure if it makes sense to have at the boolean level

  TermVec red_cube_lits, rem_cube_lits;
  reducer_.reduce_assump_unsatcore(
      formula, cube_lits, red_cube_lits, &rem_cube_lits);

  // should need some assumptions
  // formula should not be unsat on its own
  assert(red_cube_lits.size() > 0);

  pred = ic3formula_conjunction(red_cube_lits);
  // expecting a Cube here
  assert(!pred.disjunction);
}

void IC3::check_ts() const
{
  for (const auto & sv : ts_.statevars()) {
    if (sv->get_sort() != boolsort_) {
      throw PonoException("Got non-boolean state variable in bit-level IC3: "
                          + sv->to_string());
    }
  }

  for (const auto & iv : ts_.inputvars()) {
    if (iv->get_sort() != boolsort_) {
      throw PonoException("Got non-boolean input variable in bit-level IC3: "
                          + iv->to_string());
    }
  }
}

}  // namespace pono
