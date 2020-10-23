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

#include "implicit_predicate_abstractor.h"

#include "assert.h"

using namespace smt;
using namespace std;

namespace pono {

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
  Term next_pred = abs_ts_.next(pred);
  Term rel = solver_->make_term(Iff, next_pred, abstract(next_pred));

  assert(!abs_ts_.is_functional());
  abs_rts_.constrain_trans(solver_->make_term(Implies, predabs_label_, rel));
  return rel;
}

void ImplicitPredicateAbstractor::do_abstraction()
{
  // assume abs_ts_ is a perfect copy currently
  assert(abs_ts_.init() == conc_ts_.init());
  assert(abs_ts_.trans() == conc_ts_.trans());

  // create abstract variables for each next state variable
  for (auto sv : conc_ts_.statevars()) {
    Term nv = conc_ts_.next(sv);
    Term abs_nv = abs_ts_.make_statevar(nv->to_string() + "^", nv->get_sort());
    update_term_cache(nv, abs_nv);
  }

  // create the label
  assert(!predabs_label_);  // should be uninitialized
  predabs_label_ =
      abs_ts_.make_inputvar("__predabs_label", solver_->make_sort(BOOL));

  // for now assume that the abstraction is relational
  // TODO fix this -- should be able to update functional systems also
  assert(!abs_ts_.is_functional());
  // TODO: fix the population.
  // Right now state_updates, constraints, and named_terms are not updated
  Term trans = conc_ts_.trans();
  abs_rts_.set_trans(abstract(trans));
  // predabs label will force the concrete and abstract states to satisfy the
  // same predicates
  abs_rts_.constrain_inputs(predabs_label_);
}

}  // namespace pono
