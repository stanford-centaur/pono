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
  // constrain next state vars and abstract vars to agree on this predicate
  Term rel = solver_->make_term(Iff, next_pred, abstract(next_pred));
  abs_rts_.constrain_trans(rel);
  return rel;
}

void ImplicitPredicateAbstractor::do_abstraction()
{
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
}

}  // namespace pono
