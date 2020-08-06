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
    : super(ts),
      abstract_array_equality_(abstract_array_equality),
      solver_(abs_ts_.solver())
{
  do_abstraction();
}

Term ArrayAbstractor::abstract(const Term & t) const
{
  auto it = abstraction_cache_.find(t);
  if (it == abstraction_cache_.end()) {
    throw PonoException("Generic abstraction NYI");
  }
  return it->second;
}

Term ArrayAbstractor::concrete(const Term & t) const
{
  auto it = concretization_cache_.find(t);
  if (it == concretization_cache_.end()) {
    throw PonoException("Generic concretization NYI");
  }
  return it->second;
}

Term ArrayAbstractor::get_read_uf(const smt::Sort & sort) const
{
  auto it = read_ufs_.find(sort);
  if (it == read_ufs_.end()) {
    throw PonoException("No read UF found for" + sort->to_string());
  }
  return it->second;
}

Term ArrayAbstractor::get_write_uf(const smt::Sort & sort) const
{
  auto it = write_ufs_.find(sort);
  if (it == write_ufs_.end()) {
    throw PonoException("No write UF found for" + sort->to_string());
  }
  return it->second;
}

Term ArrayAbstractor::get_arrayeq_uf(const smt::Sort & sort) const
{
  // only makes sense to call if running with abstract_array_equality_
  assert(abstract_array_equality_);
  auto it = arrayeq_ufs_.find(sort);
  if (it == arrayeq_ufs_.end()) {
    throw PonoException("No write UF found for" + sort->to_string());
  }
  return it->second;
}

void ArrayAbstractor::do_abstraction() {
  abstract_vars();
  // TODO: remove these debug prints
  // cout << "S" << endl;
  // for (auto v : abs_ts_.statevars())
  // {
  //   cout << v;
  //   if (v->get_sort()->get_sort_kind() == UNINTERPRETED)
  //   {
  //     cout << ": " << get_read_uf(v->get_sort()) << endl;
  //   }
  //   else
  //   {
  //     cout << endl;
  //   }
  // }

  // cout << "I" << endl;
  // for (auto v : abs_ts_.inputvars())
  // {
  //   cout << v;
  //   if (v->get_sort()->get_sort_kind() == UNINTERPRETED)
  //   {
  //     cout << ": " << get_read_uf(v->get_sort()) << endl;
  //   }
  //   else
  //   {
  //     cout << endl;
  //   }
  // }
  // TODO: end of debug prints
  throw PonoException("got to end of implementation");
}

void ArrayAbstractor::abstract_vars()
{
  // TODO: figure out if other variables should be in cache
  //       with identity mapping (seems wasteful)
  Sort sort;
  Term abs_var;
  for (auto sv : conc_ts_.statevars()) {
    sort = sv->get_sort();
    if (sort->get_sort_kind() == ARRAY) {
      abs_var = abs_ts_.make_statevar("abs_" + sv->to_string(),
                                      abstract_array_sort(sort));
      update_term_cache(sv, abs_var);
      update_term_cache(conc_ts_.next(sv), abs_ts_.next(abs_var));
    } else {
      abs_ts_.add_statevar(sv, conc_ts_.next(sv));
    }
  }
  for (auto iv : conc_ts_.inputvars()) {
    sort = iv->get_sort();
    if (sort->get_sort_kind() == ARRAY) {
      abs_var = abs_ts_.make_statevar("abs_" + iv->to_string(),
                                      abstract_array_sort(sort));
      update_term_cache(iv, abs_var);
    } else {
      abs_ts_.add_inputvar(iv);
    }
  }
}

Sort ArrayAbstractor::abstract_array_sort(const Sort & conc_sort)
{
  assert(conc_sort->get_sort_kind() == ARRAY);
  if (abstract_sorts_.find(conc_sort) != abstract_sorts_.end()) {
    return abstract_sorts_.at(conc_sort);
  } else {
    // need to create a new uninterpreted sort
    // first: get (abstract) index and element sorts
    // TODO: figure out if getting the abstract sorts is necessary
    Sort idxsort = conc_sort->get_indexsort();
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

    // expecting UFs to not be created (this is a new abstract sort)
    assert(read_ufs_.find(abs_sort) == read_ufs_.end());
    assert(write_ufs_.find(abs_sort) == write_ufs_.end());
    assert(arrayeq_ufs_.find(abs_sort) == arrayeq_ufs_.end());

    // need to create new read and write UFs for it
    // also equality ufs if enabled
    Sort readsort =
        solver_->make_sort(FUNCTION, { abs_sort, idxsort, elemsort });
    read_ufs_[abs_sort] = solver_->make_symbol(
        "absarr_read" + std::to_string(read_ufs_.size()), readsort);

    Sort writesort =
        solver_->make_sort(FUNCTION, { abs_sort, idxsort, elemsort, abs_sort });
    write_ufs_[abs_sort] = solver_->make_symbol(
        "absarr_write" + std::to_string(write_ufs_.size()), writesort);

    if (abstract_array_equality_) {
      Sort eqsort = solver_->make_sort(
          FUNCTION, { abs_sort, abs_sort, solver_->make_sort(BOOL) });
      arrayeq_ufs_[abs_sort] = solver_->make_symbol(
          "absarr_eq" + std::to_string(arrayeq_ufs_.size()), eqsort);
    }

    return abs_sort;
  }
}

void ArrayAbstractor::update_term_cache(const Term & conc_term,
                                        const Term & abs_term)
{
  // abstraction should never change
  assert(abstraction_cache_.find(conc_term) == abstraction_cache_.end());
  assert(concretization_cache_.find(abs_term) == concretization_cache_.end());

  abstraction_cache_[conc_term] = abs_term;
  concretization_cache_[abs_term] = conc_term;
}

void ArrayAbstractor::update_sort_cache(const Sort & conc_sort,
                                        const Sort & abs_sort)
{
  // abstraction should never change
  assert(abstract_sorts_.find(conc_sort) == abstract_sorts_.end());
  assert(concrete_sorts_.find(abs_sort) == concrete_sorts_.end());

  abstract_sorts_[conc_sort] = abs_sort;
  concrete_sorts_[abs_sort] = conc_sort;
}

}  // namespace pono
