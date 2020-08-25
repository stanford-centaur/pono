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

#include "assert.h"
#include "gmpxx.h"

#include "refiners/array_axiom_enumerator.h"

using namespace smt;
using namespace std;

namespace pono {

// ArrayFinder implementation

ArrayFinder::ArrayFinder(ArrayAxiomEnumerator & aae)
    // do clear the cache -- if called again, want to add
    : aae_(aae), super(aae_.solver_, true)
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

  if (sk != ARRAY && op != Equal) {
    return Walker_Continue;
  }

  if (sk != ARRAY) {
    assert(op == Equal);
    TermVec children(term->begin(), term->end());
    assert(children.size() == 2);
    if (children[0]->get_sort()->get_sort_kind() == ARRAY) {
      Term abs_arr_eq = aae_.aa_.abstract(term);
      if (aae_.arrayeq_witnesses_.find(abs_arr_eq)
          != aae_.arrayeq_witnesses_.end()) {
        // already recorded this array equality
        return Walker_Continue;
      }

      // create a witness index for this array equality
      Term witness_idx = aae_.aa_.abs_ts().make_statevar(
          "wit_" + std::to_string(aae_.arrayeq_witnesses_.size()),
          // always uses integer index
          aae_.solver_->make_sort(INT));
      aae_.arrayeq_witnesses_[abs_arr_eq] = witness_idx;
      // add witness to index set
      aae_.index_set_.insert(witness_idx);
    }
    return Walker_Continue;
  }

  assert(sk == ARRAY);

  if (term->is_symbolic_const()) {
    // nothing to save for an array variable
    return Walker_Continue;
  }

  Term abs_term = aae_.aa_.abstract(term);
  TermVec children(term->begin(), term->end());
  TermVec abs_children(abs_term->begin(), abs_term->end());

  if (op.is_null()) {
    // constant array
    assert(children.size() == 1);
    Term val = aae_.aa_.abstract(children[0]);
    aae_.constarrs_[abs_term] = val;
  } else if (op == Store) {
    assert(abs_children.size() == 4);
    assert(children.size() == 4);
    aae_.stores_.insert(abs_term);

    // third child is index because
    // read_uf, array, index, element
    aae_.index_set_.insert(abs_children[2]);
  } else {
    assert(op == Select);
    assert(abs_children.size() == 3);
    assert(children.size() == 3);
    // third child is index, because it's
    // read_uf, array, index
    aae_.index_set_.insert(abs_children[2]);
  }

  return Walker_Continue;
}

// ArrayAxiomEnumerator implementation

ArrayAxiomEnumerator::ArrayAxiomEnumerator(Property & prop,
                                           ArrayAbstractor & aa,
                                           Unroller & un)
    : super(prop.transition_system()), prop_(prop), aa_(aa), un_(un)
{
  false_ = solver_->make_term(false);
}

bool ArrayAxiomEnumerator::enumerate_axioms(const Term & abs_trace_formula,
                                            size_t bound)
{
  // TODO: think about how to check / instantiate next-state version of axioms!!

  // clear the previous axioms
  axioms_to_check_.clear();
  violated_axioms_.clear();
  ts_axioms_.clear();
  consecutive_axioms_.clear();
  nonconsecutive_axioms_.clear();

  // Important : set bound member variable
  // used by other functions
  bound_ = bound;

  solver_->push();
  throw PonoException("NYI");
  solver_->pop();
}

// protected methods

void ArrayAxiomEnumerator::collect_arrays_and_indices()
{
  ArrayFinder af(*this);
  // just visit each (concrete) term of the transition system
  Term init = aa_.conc_ts().init();
  Term trans = aa_.conc_ts().trans();
  Term prop_term = prop_.prop();
  prop_term = aa_.concrete(prop_term);
  af.visit(init);
  af.visit(trans);
  af.visit(prop_term);

  // TODO: create lambda variables!!

  for (auto idx : index_set_) {
    if (ts_.only_curr(idx)) {
      cur_index_set_.insert(idx);
    }
  }
}

void ArrayAxiomEnumerator::check_consecutive_axioms(AxiomClass ac,
                                                    bool only_curr,
                                                    int lemma_limit)
{
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
  // TODO: fix boundary condition! if there's next in the axiom and unroll at
  // bound_
  size_t num_found_lemmas = 0;
  Term unrolled_ax;
  for (Term ax : axioms_to_check) {
    for (size_t k = 0; k <= bound_; ++k) {
      unrolled_ax = un_.at_time(ax, k);
      if (is_violated(unrolled_ax)) {
        violated_axioms_.insert(unrolled_ax);
        ts_axioms_[unrolled_ax] = ax;
        num_found_lemmas++;

        if (lemma_limit > 0 && num_found_lemmas >= lemma_limit) {
          // if given a lemma limit, then finish when that limit is reached
          return;
        }
      }
    }
  }
}

void ArrayAxiomEnumerator::check_nonconsecutive_axioms(AxiomClass ac,
                                                       bool only_curr,
                                                       int lemma_limit)
{
  throw PonoException("NYI");
}

bool ArrayAxiomEnumerator::is_violated(const Term & ax) const
{
  assert(ax->get_sort()->get_sort_kind() == BOOL);
  return solver_->get_value(ax) == false_;
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
      Sort conc_sort = aa_.concrete(idx)->get_sort();
      Term lambda = lambdas_.at(conc_sort);
      assert(lambda != idx);
      axioms_to_check.push_back(
          AxiomInstantiation(lambda_alldiff_axiom(lambda, idx), { idx }));
      if (ts_.only_curr(idx)) {
        Term next_idx = ts_.next(idx);
        assert(next_idx != lambda);
        axioms_to_check.push_back(AxiomInstantiation(
            lambda_alldiff_axiom(lambda, next_idx), { next_idx }));
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
  Term lam = lambdas_.at(conc_sort);
  Term ax = constarr_axiom(constarr, val, lam);
  assert(conc_sort->get_sort_kind() == ARRAY);
  Sort idxsort = conc_sort->get_indexsort();
  if (idxsort->get_sort_kind() == BV) {
    // IMPORTANT: if concrete index sort is finite-domain
    // need to guard with boundary conditions for soundness
    ax = solver_->make_term(Implies, lambda_guard(idxsort, lam), ax);
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
  Term lam = lambdas_.at(conc_sort);
  TermVec children(store->begin(), store->end());
  assert(children.size() == 4);  // the UF + the 3 expected arguments
  Term ax = store_read_axiom(store, lam);

  assert(conc_sort->get_sort_kind() == ARRAY);
  Sort idxsort = conc_sort->get_indexsort();
  if (idxsort->get_sort_kind() == BV) {
    // IMPORTANT: if concrete index sort is finite-domain
    // need to guard with boundary conditions for soundness
    ax = solver_->make_term(Implies, lambda_guard(idxsort, lam), ax);
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
  return solver_->make_term(Implies, arrayeq, eq_at_index);
}

Term ArrayAxiomEnumerator::arrayeq_read_lambda_axiom(const Term & arrayeq) const
{
  Term a = *(arrayeq->begin());
  Sort sort = a->get_sort();
  Sort conc_sort = aa_.concrete(sort);
  Term lam = lambdas_.at(conc_sort);
  Term ax = arrayeq_read_axiom(arrayeq, lam);
  assert(conc_sort->get_sort_kind() == ARRAY);
  Sort idxsort = conc_sort->get_indexsort();
  if (idxsort->get_sort_kind() == BV) {
    // IMPORTANT: if concrete index sort is finite-domain
    // need to guard with boundary conditions for soundness
    ax = solver_->make_term(Implies, lambda_guard(idxsort, lam), ax);
  }
  return ax;
}

Term ArrayAxiomEnumerator::lambda_alldiff_axiom(const Term & lambda,
                                                const Term & index) const
{
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

}  // namespace pono
