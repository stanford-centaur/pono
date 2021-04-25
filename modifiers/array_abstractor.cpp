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

ArrayAbstractor::AbstractionWalker::AbstractionWalker(ArrayAbstractor & aa,
                                     UnorderedTermMap * ext_cache)
    : IdentityWalker(
        aa.solver_, false, ext_cache),  // false means don't clear cache
      aa_(aa)
{
}

WalkerStepResult ArrayAbstractor::AbstractionWalker::visit_term(Term & term)
{
  if (preorder_ || in_cache(term)) {
    return Walker_Continue;
  }

  Sort sort = term->get_sort();
  SortKind sk = sort->get_sort_kind();
  Op op = term->get_op();

  TermVec cached_children;
  Term cc;
  for (auto c : term) {
    bool ok = query_cache(c, cc);
    assert(ok);  // in post-order so should always have a cache hit
    cached_children.push_back(cc);
  }

  // handle the array abstraction cases
  // constant arrays, select, store, equality
  Term res;
  if (sk == ARRAY && op.is_null()) {
    // constant array
    Term val = *(term->begin());
    string name = "constarr" + val->to_string();
    res = aa_.abs_ts_.make_statevar(name, aa_.abstract_array_sort(sort));
  } else if (op == Select) {
    assert(cached_children.size() == 2);
    Term read_uf = aa_.get_read_uf(cached_children[0]->get_sort());
    Term idx = cached_children[1];
    res = solver_->make_term(Apply, read_uf, cached_children[0], idx);
  } else if (op == Store) {
    assert(cached_children.size() == 3);
    Term write_uf = aa_.get_write_uf(aa_.abstract_array_sort(sort));
    Term idx = cached_children[1];
    res = solver_->make_term(
        Apply, { write_uf, cached_children[0], idx, cached_children[2] });
  } else if (aa_.abstract_array_equality_ && op == Equal &&
             // an equality between arrays
             (*(term->begin()))->get_sort()->get_sort_kind() == ARRAY) {
    assert(cached_children.size() == 2);
    Term arrayeq_uf = aa_.get_arrayeq_uf(cached_children[0]->get_sort());
    res = solver_->make_term(
        Apply, arrayeq_uf, cached_children[0], cached_children[1]);
  } else {
    // all arrays should be processed already except for ITE results
    // note: array variables are abstracted in a first pass earlier
    assert(sk != ARRAY || op == Ite);
    // for non-arrays, just rebuild node if necessary
    if (op.is_null()) {
      assert(!cached_children.size());
      // for values / symbols just map to itself
      res = term;
    } else {
      assert(cached_children.size());
      res = solver_->make_term(op, cached_children);
    }
  }

  assert(res);
  // use the ArrayAbstractor cache update instead of the walker's
  // it also updates the concretization cache
  aa_.update_term_cache(term, res);

  return Walker_Continue;
}

ArrayAbstractor::ConcretizationWalker::ConcretizationWalker(ArrayAbstractor & aa,
                                                            UnorderedTermMap * ext_cache)
    : IdentityWalker(
        aa.solver_, false, ext_cache),  // false means don't clear cache
      aa_(aa)
{
}

WalkerStepResult ArrayAbstractor::ConcretizationWalker::visit_term(Term & term)
{
  if (preorder_ || in_cache(term)) {
    return Walker_Continue;
  }

  Op op = term->get_op();

  // Note: all the array variables should already been in the
  // concretization cache from the initial abstraction
  // This only needs to handle the abstracted functions

  // if not an apply, then we don't need to do anything except
  // rebuild to keep changes
  if (op != Apply) {
    if (op.is_null()) {
      aa_.update_term_cache(term, term);
    } else {
      TermVec cached_children;
      Term cc;
      for (auto c : term) {
        bool ok = query_cache(c, cc);
        assert(ok);
        cached_children.push_back(cc);
      }

      Term rebuilt = solver_->make_term(op, cached_children);
      // Note: reversed order because update_term_cache takes the
      // concrete term first
      aa_.update_term_cache(rebuilt, term);
    }
    return Walker_Continue;
  }

  assert(op == Apply);
  auto it = term->begin();
  Term uf = *it;
  TermVec cached_args;
  ++it;
  while (it != term->end()) {
    Term ca;
    bool ok = query_cache(*it, ca);
    assert(ok);
    assert(ca);
    cached_args.push_back(ca);
    ++it;
  }

  Term res;
  if (aa_.read_ufs_set_.find(uf) != aa_.read_ufs_set_.end()) {
    assert(cached_args.size() == 2);
    assert(cached_args[0]->get_sort()->get_sort_kind() == ARRAY);
    Term idx = cached_args[1];
    res = solver_->make_term(Select, cached_args[0], idx);
  } else if (aa_.write_ufs_set_.find(uf) != aa_.write_ufs_set_.end()) {
    assert(cached_args.size() == 3);
    assert(cached_args[0]->get_sort()->get_sort_kind() == ARRAY);
    Term idx = cached_args[1];
    res = solver_->make_term(Store, cached_args[0], idx, cached_args[1]);
  } else if (aa_.arrayeq_ufs_set_.find(uf) != aa_.arrayeq_ufs_set_.end()) {
    assert(cached_args.size() == 2);
    assert(cached_args[0]->get_sort()->get_sort_kind() == ARRAY);
    assert(cached_args[1]->get_sort()->get_sort_kind() == ARRAY);
    res = solver_->make_term(Equal, cached_args);
  } else {
    // just rebuild. this is not an abstract uf
    cached_args.insert(cached_args.begin(), uf);
    res = solver_->make_term(op, cached_args);
  }

  assert(res);
  // Note: reversed order because update_term_cache takes the
  // concrete term first
  aa_.update_term_cache(res, term);

  return Walker_Continue;
}

ArrayAbstractor::ArrayAbstractor(const TransitionSystem & conc_ts,
                                 TransitionSystem & abs_ts,
                                 bool abstract_array_equality)
    : super(conc_ts, abs_ts),
      abstract_array_equality_(abstract_array_equality),
      solver_(abs_ts_.solver()),
      abs_walker_(*this, &abstraction_cache_),
      conc_walker_(*this, &concretization_cache_)
{
}

Term ArrayAbstractor::abstract(Term & t) { return abs_walker_.visit(t); }

Term ArrayAbstractor::concrete(Term & t) { return conc_walker_.visit(t); }

Sort ArrayAbstractor::abstract(Sort & s)
{
  auto it = abstract_sorts_.find(s);
  if (it != abstract_sorts_.end()) {
    return it->second;
  } else {
    return s;
  }
}

Sort ArrayAbstractor::concrete(Sort & s)
{
  auto it = concrete_sorts_.find(s);
  if (it != concrete_sorts_.end()) {
    return it->second;
  } else {
    return s;
  }
}

Term ArrayAbstractor::get_read_uf(const smt::Sort & sort) const
{
  auto it = read_ufs_.find(sort);
  if (it == read_ufs_.end()) {
    throw PonoException("No read UF found for " + sort->to_string());
  }
  return it->second;
}

Term ArrayAbstractor::get_write_uf(const smt::Sort & sort) const
{
  auto it = write_ufs_.find(sort);
  if (it == write_ufs_.end()) {
    throw PonoException("No write UF found for " + sort->to_string());
  }
  return it->second;
}

Term ArrayAbstractor::get_arrayeq_uf(const smt::Sort & sort) const
{
  // only makes sense to call if running with abstract_array_equality_
  assert(abstract_array_equality_);
  auto it = arrayeq_ufs_.find(sort);
  if (it == arrayeq_ufs_.end()) {
    throw PonoException("No array equality abstraction found for: "
                        + sort->to_string());
  }
  return it->second;
}

void ArrayAbstractor::do_abstraction()
{
  abstract_vars();
  Term init = conc_ts_.init();
  Term trans = conc_ts_.trans();
  Term abs_init = abstract(init);
  Term abs_trans = abstract(trans);

  // need a relational system
  // but generic abstractor does not require a relational system
  // do a cast
  assert(!abs_ts_.is_functional());
  RelationalTransitionSystem & abs_rts =
      static_cast<RelationalTransitionSystem &>(abs_ts_);
  // the calls to abstract have already added the
  // (possibly abstracted) variables
  // now we just need to set the initial states and trans
  abs_rts.set_init(abs_init);
  abs_rts.set_trans(abs_trans);
  // TODO: constraints and state_updates are not updated correctly
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
      abs_var = abs_ts_.make_inputvar("abs_" + iv->to_string(),
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
    Term read_uf = solver_->make_symbol(
        "absarr_read" + std::to_string(read_ufs_.size()), readsort);
    read_ufs_[abs_sort] = read_uf;
    read_ufs_set_.insert(read_uf);

    Sort writesort =
        solver_->make_sort(FUNCTION, { abs_sort, idxsort, elemsort, abs_sort });
    Term write_uf = solver_->make_symbol(
        "absarr_write" + std::to_string(write_ufs_.size()), writesort);

    write_ufs_[abs_sort] = write_uf;
    write_ufs_set_.insert(write_uf);

    if (abstract_array_equality_) {
      Sort eqsort = solver_->make_sort(
          FUNCTION, { abs_sort, abs_sort, solver_->make_sort(BOOL) });
      Term arrayeq = solver_->make_symbol(
          "absarr_eq" + std::to_string(arrayeq_ufs_.size()), eqsort);
      arrayeq_ufs_[abs_sort] = arrayeq;
      arrayeq_ufs_set_.insert(arrayeq);
    }

    return abs_sort;
  }
}

void ArrayAbstractor::update_sort_cache(const Sort & conc_sort,
                                        const Sort & abs_sort)
{
  // sort abstraction should never change
  assert(abstract_sorts_.find(conc_sort) == abstract_sorts_.end());
  assert(concrete_sorts_.find(abs_sort) == concrete_sorts_.end());

  abstract_sorts_[conc_sort] = abs_sort;
  concrete_sorts_[abs_sort] = conc_sort;
}

}  // namespace pono
