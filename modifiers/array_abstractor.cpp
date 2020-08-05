/*********************                                                  */
/*! \file array_abstractor.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract arrays using uninterpreted functions.
**
**
**/

#include "assert.h"

#include "modifiers/array_abstractor.h"
#include "utils/exceptions.h"

using namespace smt;
using namespace std;

namespace pono {

ArrayAbstractor::ArrayAbstractor(const TransitionSystem & ts,
                                 bool abstract_array_equality)
    : super(ts), abstract_array_equality_(abstract_array_equality)
{
  do_abstraction();
}

Term ArrayAbstractor::abstract(const Term & t) const
{
  throw PonoException("NYI");
}

Term ArrayAbstractor::concrete(const Term & t) const
{
  throw PonoException("NYI");
}

void ArrayAbstractor::do_abstraction() { throw PonoException("NYI"); }

void ArrayAbstractor::abstract_array_vars()
{
  Sort sort;
  Term abs_var;
  for (auto sv : conc_ts_.statevars()) {
    sort = sv->get_sort();
    if (sort->get_sort_kind() == ARRAY) {
      abs_var = abs_ts_->make_statevar("abs_" + sv->to_string(),
                                       abstract_array_sort(sort));
      update_term_cache(sv, abs_var);
    }
  }
  for (auto iv : conc_ts_.inputvars()) {
    sort = iv->get_sort();
    if (sort->get_sort_kind() == ARRAY) {
      abs_var = abs_ts_->make_statevar("abs_" + iv->to_string(),
                                       abstract_array_sort(sort));
      update_term_cache(iv, abs_var);
    }
  }
}

Sort abstract_array_sort(const Sort & conc_sort)
{
  assert(conc_sort->get_sort_kind() == ARRAY);
  if (abstract_sorts_.find(conc_sort) != abstract_sorts_.end()) {
    return abstract_sorts_.at(conc_sort);
  } else {
    // need to create a new uninterpreted sort
    // first: get (abstract) index and element sorts
    Sort idxsort = conc_sort->get_indexsort();
    // TODO: figure out if getting the abstract sorts is necessary
    if (idxsort->get_sort_kind() == ARRAY) {
      idxsort = abstract_array_sort(idxsort);
    }
    Sort elemsort = conc_sort->get_elemsort();
    if (elemsort->get_sort_kind() == ARRAY) {
      elemsort = abstract_array_sort(elemsort);
    }

    Sort abs_sort = solver_->make_sort(
        "Array_" + idxsort->to_string() + "_" + elemsort->to_string(), 0);
    update_sort_cache(conc_sort, abs_sort);
  }
}

void update_term_cache(const Term & conc_term, const Term & abs_term)
{
  // abstraction should never change
  assert(abstraction_cache_.find(conc_term) == abstraction_cache_.end());
  assert(concretization_cache_.find(abs_term) == concretization_cache_.end());

  abstraction_cache_[conc_term] = abs_term;
  concretization_cache_[abs_term] = conc_term;
}

void update_sort_cache(const Sort & conc_sort, const Sort & abs_sort)
{
  // abstraction should never change
  assert(abstract_sorts_.find(conc_sort) == abstract_sorts_.end());
  assert(concrete_sorts_.find(abs_sort) == concrete_sorts_.end());

  abstract_sorts_[conc_sort] = abs_sort;
  concrete_sorts_[abs_sort] = conc_sort;
}

}  // namespace pono
