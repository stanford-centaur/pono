/*********************                                                  */
/*! \file array_axiom_enumerator.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Class for enumerating array axioms over an array abstraction
**        produced by ArrayAbstractor (see array_abstractor.[h, cpp])
**
**/

#include "refiners/array_axiom_enumerator.h"

#include <cassert>

#include "gmpxx.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

string to_string(AxiomClass ac)
{
  switch (ac) {
    case CONSTARR: return "CONSTARR";
    case CONSTARR_LAMBDA: return "CONSTARR_LAMBDA";
    case STORE_WRITE: return "STORE_WRITE";
    case STORE_READ: return "STORE_READ";
    case STORE_READ_LAMBDA: return "STORE_READ_LAMBDA";
    case ARRAYEQ_WITNESS: return "ARRAYEQ_WITNESS";
    case ARRAYEQ_READ: return "ARRAYEQ_READ";
    case ARRAYEQ_READ_LAMBDA: return "ARRAYEQ_READ_LAMBDA";
    case LAMBDA_ALLDIFF: return "LAMBDA_ALLDIFF";
    default: throw PonoException("Unhandled AxiomClass in to_string");
  }
}

// ArrayFinder implementation

ArrayFinder::ArrayFinder(ArrayAxiomEnumerator & aae)
    // do clear the cache -- if called again, want to add
    : aae_(aae), super(aae.solver_, true)
{
}

WalkerStepResult ArrayFinder::visit_term(Term & term)
{
  if (!preorder_) {
    return Walker_Continue;
  }

  // just an identity mapping
  // needed for correct term traversal
  save_in_cache(term, term);

  Sort sort = term->get_sort();
  SortKind sk = sort->get_sort_kind();
  Op op = term->get_op();

  if ((sk != ARRAY && op != Select && op != Equal) || op == Ite) {
    // if not an array or array equality, then nothing to save
    // similarly, don't need to save anything for an Ite,
    // even if it is an array sort
    return Walker_Continue;
  }

  Term abs_term = aae_.aa_.abstract(term);
  TermVec children(term->begin(), term->end());
  TermVec abs_children(abs_term->begin(), abs_term->end());

  if (sk != ARRAY && op == Equal) {
    assert(children.size() == 2);
    if (children[0]->get_sort()->get_sort_kind() == ARRAY) {
      Term abs_arr_eq = aae_.aa_.abstract(term);
      if (aae_.arrayeq_witnesses_.find(abs_arr_eq)
          != aae_.arrayeq_witnesses_.end()) {
        // already recorded this array equality
        return Walker_Continue;
      }

      Sort idxsort = children[0]->get_sort()->get_indexsort();

      // create a witness index for this array equality
      Term witness_idx = aae_.aa_.abs_ts().make_statevar(
          "wit_" + std::to_string(aae_.arrayeq_witnesses_.size()), idxsort);
      aae_.arrayeq_witnesses_[abs_arr_eq] = witness_idx;
      aae_.witnesses_to_idxsort_[witness_idx] =
          children[0]->get_sort()->get_indexsort();
      // add witness to index set
      aae_.index_set_.insert(witness_idx);
    }
    return Walker_Continue;
  } else if (sk != ARRAY && op == Select) {
    assert(op == Select);
    assert(abs_children.size() == 3);
    assert(children.size() == 2);
    // third child is index, because it's
    // read_uf, array, index
    aae_.index_set_.insert(abs_children[2]);
    return Walker_Continue;
  }

  assert(sk == ARRAY);

  if (term->is_symbolic_const()) {
    // nothing to save for an array variable
    return Walker_Continue;
  } else if (op.is_null()) {
    // constant array
    assert(children.size() == 1);
    Term val = aae_.aa_.abstract(children[0]);
    aae_.constarrs_[abs_term] = val;
  } else if (op == Store) {
    assert(abs_children.size() == 4);
    assert(children.size() == 3);
    aae_.stores_.insert(abs_term);

    // third child is index because
    // read_uf, array, index, element
    aae_.index_set_.insert(abs_children[2]);
  } else {
    // this should be all the cases
    assert(false);
  }

  return Walker_Continue;
}

// ArrayAxiomEnumerator implementation

ArrayAxiomEnumerator::ArrayAxiomEnumerator(ArrayAbstractor & aa,
                                           Unroller & un,
                                           const Term & prop,
                                           bool red_axioms)
    : super(aa.abs_ts()), aa_(aa), un_(un), reduce_axioms_unsatcore_(red_axioms)
{
  conc_bad_ = solver_->make_term(Not, prop);
  false_ = solver_->make_term(false);
}

void ArrayAxiomEnumerator::initialize()
{
  assert(!initialized_);

  initialized_ = true;

  collect_arrays_and_indices();
  create_lambda_indices();
}

bool ArrayAxiomEnumerator::enumerate_axioms(const Term & abs_trace_formula,
                                            size_t bound,
                                            bool include_nonconsecutive)
{
  return enumerate_axioms(
      abs_trace_formula, bound, include_nonconsecutive, false);
}

bool ArrayAxiomEnumerator::enumerate_axioms(const Term & abs_trace_formula,
                                            size_t bound,
                                            bool include_nonconsecutive,
                                            bool skip_lambda_axioms)
{
  assert(initialized_);

  // IMPORTANT: clear state from last run
  clear_state();

  // TODO: think about how to check / instantiate next-state version of axioms!!
  // Important : set bound member variable
  // used by other functions
  bound_ = bound;

  // don't want to pollute state of solver
  // do all solving in a new context
  solver_->push();

  solver_->assert_formula(abs_trace_formula);
  Result res = solver_->check_sat();
  TermVec all_violated_axioms;

  // use only current axioms if the bound is zero
  // e.g. only the initial state
  // this is because we can only add axioms over current state variables
  // to initial states
  bool only_curr = (bound == 0);
  while (res.is_sat()) {
    bool found_lemmas = false;

    // check axioms
    // heuristic order -- all need to be checked for completeness
    // but might not need to add all of them to prove a property
    // preferring axioms that don't enumerate indices first
    // except not lambda axioms -- those are fairly rare

    // experimenting with new heuristic order
    found_lemmas |= check_consecutive_axioms(ARRAYEQ_READ, only_curr);
    if (!skip_lambda_axioms) {
      found_lemmas |= check_consecutive_axioms(ARRAYEQ_READ_LAMBDA, only_curr);
    }
    found_lemmas |= check_consecutive_axioms(STORE_WRITE, only_curr);
    found_lemmas |= check_consecutive_axioms(STORE_READ, only_curr);
    if (!skip_lambda_axioms) {
      found_lemmas |= check_consecutive_axioms(STORE_READ_LAMBDA, only_curr);
    }
    found_lemmas |= check_consecutive_axioms(CONSTARR, only_curr);
    if (!skip_lambda_axioms) {
      found_lemmas |= check_consecutive_axioms(CONSTARR_LAMBDA, only_curr);
    }
    found_lemmas |= check_consecutive_axioms(ARRAYEQ_WITNESS, only_curr);

    if (!found_lemmas && !skip_lambda_axioms) {
      // NOTE: don't need non-consecutive version of these axioms
      //       all different over current and next is sufficient to be all
      //       different for all time
      found_lemmas |= check_consecutive_axioms(LAMBDA_ALLDIFF, only_curr);
    }

    // END EXP

    // check non-consecutive axioms now if no other lemmas have been found
    // need to check at unrolled indices
    // for performance, we prefer indices that are a short distance
    // from the property violation (at bound_)
    // this will result in less auxiliary variables to make the axiom
    // consecutive
    int k = bound_;
    while (include_nonconsecutive && !found_lemmas && k >= 0) {
      found_lemmas |= check_nonconsecutive_axioms(ARRAYEQ_READ, only_curr, k);
      found_lemmas |= check_nonconsecutive_axioms(STORE_READ, only_curr, k);
      found_lemmas |= check_nonconsecutive_axioms(CONSTARR, only_curr, k);
      k--;
    }

    if (!found_lemmas) {
      // there appears to be a concrete counterexample
      solver_->pop();
      return false;
    }

    assert(violated_axioms_.size());
    for (auto ax : violated_axioms_) {
      // save the axiom
      all_violated_axioms.push_back(ax);
      if (!reduce_axioms_unsatcore_) {
        solver_->assert_formula(ax);
      }
    }
    // reset violated_axioms_ so we don't add the same axioms again
    violated_axioms_.clear();

    if (reduce_axioms_unsatcore_) {
      res = solver_->check_sat_assuming(all_violated_axioms);
    } else {
      res = solver_->check_sat();
    }
    assert(!res.is_unknown());  // this algorithm assumes decidable theories
  }

  assert(res.is_unsat());  // ruled out the trace

  UnorderedTermSet core_set;
  if (reduce_axioms_unsatcore_) {
    try {
      solver_->get_unsat_assumptions(core_set);
      logger.log(1,
                 "Reduced bmc axioms: {}/{}",
                 core_set.size(),
                 all_violated_axioms.size());
    }
    catch (SmtException & e) {
      // if core is empty, that's fine -- continue
      // but it should only happen if no axioms were added
      assert(!all_violated_axioms.size());
    }
  }

  // populate axioms
  for (auto ax : all_violated_axioms) {
    if (reduce_axioms_unsatcore_ && core_set.find(ax) == core_set.end()) {
      // if not in unsat core
      // then don't keep this axiom
      continue;
    }
    if (ts_axioms_.find(ax) != ts_axioms_.end()) {
      // this is a consecutive axiom
      consecutive_axioms_.insert(ts_axioms_.at(ax));
    } else {
      // this is a non-consecutive axiom
      // expect it in the map to AxiomInstantiations
      assert(to_axiom_inst_.find(ax) != to_axiom_inst_.end());
      nonconsecutive_axioms_.push_back(to_axiom_inst_.at(ax));
    }
  }

  solver_->pop();
  return true;
}

void ArrayAxiomEnumerator::add_index(const Term & idx)
{
  index_set_.insert(idx);
  if (aa_.abs_ts().only_curr(idx)) {
    cur_index_set_.insert(idx);
  }
}

bool ArrayAxiomEnumerator::check_consecutive_axioms(AxiomClass ac,
                                                    bool only_curr,
                                                    int lemma_limit)
{
  assert(initialized_);

  logger.log(3, "Checking consecutive axioms for class: {}", to_string(ac));
  UnorderedTermSet & indices = only_curr ? cur_index_set_ : index_set_;

  UnorderedTermSet axioms_to_check;
  if (index_axiom_classes.find(ac) == index_axiom_classes.end()) {
    axioms_to_check = non_index_axioms(ac);
  } else {
    for (AxiomInstantiation ax_inst : index_axioms(ac, indices)) {
      // consecutive axioms don't need to keep track of specific index
      // instantiations just keep the axiom term
      axioms_to_check.insert(ax_inst.ax);
    }
  }

  // TODO: make sure we're covering the current/next
  //       version of axioms correctly
  //       not explicitly calling next here -- is that a problem?
  //       should be okay as long as we add next version of axioms
  //       that are only over state variables
  size_t num_found_lemmas = 0;
  Term unrolled_ax;
  for (Term ax : axioms_to_check) {
    if (only_curr && !ts_.only_curr(ax)) {
      // if requesting axioms over only current state variables
      // and this axiom isn't, then continue
      continue;
    }
    // bound to check until depends on whether there are inputs/next state vars
    // in the axiom
    size_t max_k = ts_.only_curr(ax) ? bound_ : bound_ - 1;
    for (size_t k = 0; k <= max_k; ++k) {
      unrolled_ax = un_.at_time(ax, k);
      if (is_violated(unrolled_ax)) {
        logger.log(4, "Violated Axiom: {}", unrolled_ax);
        violated_axioms_.insert(unrolled_ax);
        ts_axioms_[unrolled_ax] = ax;
        num_found_lemmas++;

        if (lemma_limit > 0 && num_found_lemmas >= lemma_limit) {
          // if given a lemma limit, then finish when that limit is reached
          return num_found_lemmas;
        }
      }
    }
  }

  return num_found_lemmas;
}

bool ArrayAxiomEnumerator::check_nonconsecutive_axioms(AxiomClass ac,
                                                       bool only_curr,
                                                       size_t i,
                                                       int lemma_limit)
{
  assert(initialized_);

  logger.log(3,
             "Checking nonconsecutive axioms at {} for class: {}",
             i,
             to_string(ac));
  // there are no non-consecutive axioms that don't instantiate axioms
  // thus the AxiomClass must be one parameterized by an index
  assert(index_axiom_classes.find(ac) != index_axiom_classes.end());

  // must be within bound
  assert(i <= bound_);

  UnorderedTermSet & indices = only_curr ? cur_index_set_ : index_set_;
  UnorderedTermSet unrolled_indices;
  Term unrolled_idx;
  for (auto idx : indices) {
    if (i == bound_ && !ts_.only_curr(idx)) {
      // IMPORTANT: cannot instantiate anything containing
      // an input variable at the bound -- because then
      // the target would not be a state variable
      // since there is no delay
      // this is consistent with transition system formulation
      // it doesn't make sense to talk about the next of an input
      // and the bound is essentially the last next unrolling
      continue;
    }
    unrolled_idx = un_.at_time(idx, i);
    unrolled_indices.insert(unrolled_idx);
    untime_index_cache_[unrolled_idx] = idx;
  }

  // check these axioms
  // Note: using staged unrolling -- i.e. indices already unrolled
  // but the rest of the axiom is not, until later
  size_t num_found_lemmas = 0;
  Term unrolled_ax;
  for (AxiomInstantiation ax_inst : index_axioms(ac, unrolled_indices)) {
    if (only_curr && !ts_.only_curr(ax_inst.ax)) {
      // if requesting axioms over only current state variables
      // and this axiom isn't, then continue
      continue;
    }
    // bound to check until depends on whether there are inputs/next state vars
    // in the axiom
    // Note: the ax_inst already has an unrolled index, but the other variables
    // have not been unrolled. Example,
    //    ax_inst: i@1 != j -> read(write(a, j, e), i@1) = read(a, i@1)
    //    this would get unrolled for different k values, for k = 3 it would be:
    //             i@1 != j@3 -> read(write(a@3, j@3, e@3), i@1) = read(a@3,
    //             i@1)
    size_t max_k = ts_.only_curr(ax_inst.ax) ? bound_ : bound_ - 1;
    for (size_t k = 0; k <= max_k; ++k) {
      unrolled_ax = un_.at_time(ax_inst.ax, k);
      if (is_violated(unrolled_ax)) {
        logger.log(4, "Violated NonConsecutive Axiom: {}", unrolled_ax);
        violated_axioms_.insert(unrolled_ax);
        to_axiom_inst_.insert({ unrolled_ax, ax_inst });
        num_found_lemmas++;

        if (lemma_limit > 0 && num_found_lemmas >= lemma_limit) {
          // if given a lemma limit, then finish when that limit is reached
          return num_found_lemmas;
        }
      }
    }
  }
  return num_found_lemmas;
}

// protected methods

void ArrayAxiomEnumerator::collect_arrays_and_indices()
{
  assert(initialized_);

  ArrayFinder af(*this);
  // just visit each (concrete) term of the transition system
  Term init = aa_.conc_ts().init();
  Term trans = aa_.conc_ts().trans();
  af.visit(init);
  af.visit(trans);
  assert(conc_bad_);
  af.visit(conc_bad_);

  for (auto idx : index_set_) {
    if (ts_.only_curr(idx)) {
      cur_index_set_.insert(idx);
    }
  }
}

void ArrayAxiomEnumerator::create_lambda_indices()
{
  unordered_set<Sort> conc_array_idx_sorts;
  const TransitionSystem & conc_ts = aa_.conc_ts();
  for (auto sv : conc_ts.statevars()) {
    Sort sort = sv->get_sort();
    if (sort->get_sort_kind() == ARRAY) {
      conc_array_idx_sorts.insert(sort->get_indexsort());
    }
  }
  for (auto i : conc_ts.inputvars()) {
    Sort sort = i->get_sort();
    if (sort->get_sort_kind() == ARRAY) {
      conc_array_idx_sorts.insert(sort->get_indexsort());
    }
  }

  // create a lambda var for each sort
  Sort intsort = solver_->make_sort(INT);
  size_t lam_num = 0;
  // TODO: make this less confusing
  // ts_ is const
  // but the same TransitionSystem is accessible through the ArrayAbstractor
  // and is mutable
  TransitionSystem & mutable_ts = aa_.abs_ts();
  for (auto idxsort : conc_array_idx_sorts) {
    Term lam = mutable_ts.make_statevar("lambda_" + std::to_string(lam_num++),
                                        // always using an integer sort for
                                        // lambdas to avoid finite domain issues
                                        intsort);
    lambdas_[idxsort] = lam;
  }
}

void ArrayAxiomEnumerator::clear_state()
{
  // clear the previous axioms
  axioms_to_check_.clear();
  violated_axioms_.clear();
  ts_axioms_.clear();
  to_axiom_inst_.clear();
  consecutive_axioms_.clear();
  nonconsecutive_axioms_.clear();
}

bool ArrayAxiomEnumerator::is_violated(const Term & ax) const
{
  assert(ax->get_sort()->get_sort_kind() == BOOL);
  Term val = solver_->get_value(ax);
  // expecting to have solver options such that
  // we can get the value of arbitrary axioms
  assert(val->is_value());
  return val == false_;
}

UnorderedTermSet ArrayAxiomEnumerator::non_index_axioms(AxiomClass ac)
{
  assert(index_axiom_classes.find(ac) == index_axiom_classes.end());

  UnorderedTermSet axioms_to_check;
  if (ac == CONSTARR_LAMBDA) {
    for (auto elem : constarrs_) {
      axioms_to_check.insert(constarr_lambda_axiom(elem.first, elem.second));
    }
  } else if (ac == STORE_WRITE) {
    for (auto st : stores_) {
      axioms_to_check.insert(store_write_axiom(st));
    }
  } else if (ac == STORE_READ_LAMBDA) {
    for (auto st : stores_) {
      axioms_to_check.insert(store_read_lambda_axiom(st));
    }
  } else if (ac == ARRAYEQ_WITNESS) {
    for (auto elem : arrayeq_witnesses_) {
      axioms_to_check.insert(arrayeq_witness_axiom(elem.first));
    }
  } else if (ac == ARRAYEQ_READ_LAMBDA) {
    for (auto elem : arrayeq_witnesses_) {
      axioms_to_check.insert(arrayeq_read_lambda_axiom(elem.first));
    }
  } else {
    throw PonoException("Unhandled AxiomClass");
  }

  return axioms_to_check;
}

AxiomVec ArrayAxiomEnumerator::index_axioms(AxiomClass ac,
                                            UnorderedTermSet & indices)
{
  assert(index_axiom_classes.find(ac) != index_axiom_classes.end());

  AxiomVec axioms_to_check;
  for (auto idx : indices) {
    if (ac == CONSTARR) {
      for (auto elem : constarrs_) {
        axioms_to_check.push_back(AxiomInstantiation(
            constarr_axiom(elem.first, elem.second, idx), { idx }));
      }
    } else if (ac == STORE_READ) {
      for (auto st : stores_) {
        axioms_to_check.push_back(
            AxiomInstantiation(store_read_axiom(st, idx), { idx }));
      }
    } else if (ac == ARRAYEQ_READ) {
      for (auto elem : arrayeq_witnesses_) {
        axioms_to_check.push_back(
            AxiomInstantiation(arrayeq_read_axiom(elem.first, idx), { idx }));
      }
    } else if (ac == LAMBDA_ALLDIFF) {
      // CRUCIAL: only instantiate all different axioms over matching sorts
      // e.g. the lambda instantiated for a particular index sort
      // can look up the lambda by the (concrete) index sort
      Term lam = get_lambda(idx);
      assert(lam != idx);
      axioms_to_check.push_back(
          AxiomInstantiation(lambda_alldiff_axiom(lam, idx), { idx }));
      if (ts_.only_curr(idx)) {
        Term next_idx = ts_.next(idx);
        assert(next_idx != lam);
        axioms_to_check.push_back(AxiomInstantiation(
            lambda_alldiff_axiom(lam, next_idx), { next_idx }));
      }
    } else {
      throw PonoException("Unhandled AxiomClass");
    }
  }
  return axioms_to_check;
}

Term ArrayAxiomEnumerator::constarr_axiom(const Term & constarr,
                                          const Term & val,
                                          const Term & index) const
{
  Term read_uf = aa_.get_read_uf(constarr->get_sort());
  return solver_->make_term(
      Equal, solver_->make_term(Apply, read_uf, constarr, index), val);
}

Term ArrayAxiomEnumerator::constarr_lambda_axiom(const Term & constarr,
                                                 const Term & val) const
{
  Sort sort = constarr->get_sort();
  Term read_uf = aa_.get_read_uf(sort);
  Sort conc_sort = aa_.concrete(sort);
  assert(conc_sort->get_sort_kind() == ARRAY);
  Sort conc_idx_sort = conc_sort->get_indexsort();
  Term lam = lambdas_.at(conc_idx_sort);
  Term casted_lam = cast_lambda(conc_idx_sort, lam);
  Term ax = constarr_axiom(constarr, val, casted_lam);
  if (conc_idx_sort->get_sort_kind() == BV) {
    // IMPORTANT: if concrete index sort is finite-domain
    // need to guard with boundary conditions for soundness
    ax = solver_->make_term(Implies, lambda_guard(conc_idx_sort, lam), ax);
  }
  return ax;
}

Term ArrayAxiomEnumerator::store_write_axiom(const Term & store) const
{
  Term read_uf = aa_.get_read_uf(store->get_sort());
  TermVec children(store->begin(), store->end());
  assert(children.size() == 4);  // the UF + the 3 expected arguments
  Term idx = children[2];
  Term val = children[3];
  return solver_->make_term(
      Equal, solver_->make_term(Apply, read_uf, store, idx), val);
}

Term ArrayAxiomEnumerator::store_read_axiom(const Term & store,
                                            const Term & index) const
{
  Term read_uf = aa_.get_read_uf(store->get_sort());
  TermVec children(store->begin(), store->end());
  assert(children.size() == 4);  // the UF + the 3 expected arguments
  Term a = children[1];
  Term wr_idx = children[2];
  Term val = children[3];
  Term antecedent = solver_->make_term(Distinct, index, wr_idx);
  Term read_store = solver_->make_term(Apply, read_uf, store, index);
  Term read_a = solver_->make_term(Apply, read_uf, a, index);
  return solver_->make_term(
      Implies, antecedent, solver_->make_term(Equal, read_store, read_a));
}

Term ArrayAxiomEnumerator::store_read_lambda_axiom(const Term & store) const
{
  Sort sort = store->get_sort();
  Term read_uf = aa_.get_read_uf(sort);
  Sort conc_sort = aa_.concrete(sort);
  assert(conc_sort->get_sort_kind() == ARRAY);
  Sort conc_idx_sort = conc_sort->get_indexsort();
  Term lam = lambdas_.at(conc_idx_sort);
  Term casted_lam = cast_lambda(conc_idx_sort, lam);
  TermVec children(store->begin(), store->end());
  assert(children.size() == 4);  // the UF + the 3 expected arguments
  Term ax = store_read_axiom(store, casted_lam);

  if (conc_idx_sort->get_sort_kind() == BV) {
    // IMPORTANT: if concrete index sort is finite-domain
    // need to guard with boundary conditions for soundness
    ax = solver_->make_term(Implies, lambda_guard(conc_idx_sort, lam), ax);
  }
  return ax;
}

Term ArrayAxiomEnumerator::arrayeq_witness_axiom(const Term & arrayeq) const
{
  Term witness = arrayeq_witnesses_.at(arrayeq);
  TermVec children(arrayeq->begin(), arrayeq->end());
  Term a, b;
  if (aa_.abstract_array_equality()) {
    assert(children.size() == 3);  // the UF + 2 arrays
    assert(arrayeq->get_op() == Apply);
    a = children[1];
    b = children[2];
  } else {
    assert(children.size() == 2);
    assert(arrayeq->get_op() == Equal);
    a = children[0];
    b = children[1];
  }
  assert(a);
  assert(b);

  Term read_uf = aa_.get_read_uf(a->get_sort());
  Term eq_at_witness =
      solver_->make_term(Equal,
                         solver_->make_term(Apply, read_uf, a, witness),
                         solver_->make_term(Apply, read_uf, b, witness));
  return solver_->make_term(Implies, eq_at_witness, arrayeq);
}

Term ArrayAxiomEnumerator::arrayeq_read_axiom(const Term & arrayeq,
                                              const Term & index) const
{
  TermVec children(arrayeq->begin(), arrayeq->end());
  Term a, b;
  if (aa_.abstract_array_equality()) {
    assert(children.size() == 3);  // the UF + 2 arrays
    assert(arrayeq->get_op() == Apply);
    a = children[1];
    b = children[2];
  } else {
    assert(children.size() == 2);
    assert(arrayeq->get_op() == Equal);
    a = children[0];
    b = children[1];
  }
  assert(a);
  assert(b);

  Term read_uf = aa_.get_read_uf(a->get_sort());
  Term eq_at_index =
      solver_->make_term(Equal,
                         solver_->make_term(Apply, read_uf, a, index),
                         solver_->make_term(Apply, read_uf, b, index));
  return solver_->make_term(Implies,
                            solver_->make_term(Not, eq_at_index),
                            solver_->make_term(Not, arrayeq));
}

Term ArrayAxiomEnumerator::arrayeq_read_lambda_axiom(const Term & arrayeq) const
{
  // TODO: fix the const issues for passing to abstractor
  Term ae = arrayeq;
  Term conc_eq = aa_.concrete(ae);
  Term conc_a = *(conc_eq->begin());
  Sort conc_sort = conc_a->get_sort();
  assert(conc_sort->get_sort_kind() == ARRAY);
  Sort conc_idx_sort = conc_sort->get_indexsort();
  Term lam = lambdas_.at(conc_idx_sort);
  Term casted_lam = cast_lambda(conc_idx_sort, lam);
  Term ax = arrayeq_read_axiom(arrayeq, casted_lam);
  if (conc_idx_sort->get_sort_kind() == BV) {
    // IMPORTANT: if concrete index sort is finite-domain
    // need to guard with boundary conditions for soundness
    ax = solver_->make_term(Implies, lambda_guard(conc_idx_sort, lam), ax);
  }
  return ax;
}

Term ArrayAxiomEnumerator::lambda_alldiff_axiom(const Term & lambda,
                                                const Term & index) const
{
  if (lambda->get_sort() != index->get_sort()) {
    if (index->get_sort()->get_sort_kind() == BV) {
      Term int_index = solver_->make_term(BV_To_Nat, index);
      return solver_->make_term(Distinct, lambda, int_index);
    } else {
      throw PonoException("Unsupported index support for lambda comparison");
    }
  }
  return solver_->make_term(Distinct, lambda, index);
}

Term ArrayAxiomEnumerator::lambda_guard(const Sort & sort,
                                        const Term & lam) const
{
  Sort intsort = solver_->make_sort(INT);

  // wrote this with bit-vectors in mind
  // won't work out-of-the-box for other sorts
  assert(sort->get_sort_kind() == BV);
  // lambda should be abstracted as an integer
  assert(lam->get_sort() == intsort);

  Term zero = solver_->make_term(0, intsort);
  mpz_class maxval_gmp(std::string(sort->get_width(), '1'), 2);
  Term maxval = solver_->make_term(maxval_gmp.get_str(10), intsort);

  return solver_->make_term(And,
                            solver_->make_term(Ge, lam, zero),
                            solver_->make_term(Le, lam, maxval));
}

Term ArrayAxiomEnumerator::get_lambda(Term idx)
{
  if (witnesses_to_idxsort_.find(idx) != witnesses_to_idxsort_.end()) {
    // witness index does not have a concrete version
    // (only appears in the abstract system)
    return lambdas_.at(witnesses_to_idxsort_.at(idx));
  } else {
    Term conc_idx = aa_.concrete(idx);
    Sort conc_sort = aa_.concrete(idx)->get_sort();
    return lambdas_.at(conc_sort);
  }
}

Term ArrayAxiomEnumerator::cast_lambda(const Sort & sort,
                                       const Term & lam) const
{
  if (lam->get_sort() != sort) {
    if (sort->get_sort_kind() == BV) {
      return solver_->make_term(Op(Int_To_BV, sort->get_width()), lam);
    } else {
      throw PonoException("Unhandled sort in cast_lambda");
    }
  }
  return lam;
}

}  // namespace pono
