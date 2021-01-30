/*********************                                                        */
/*! \file term_walkers.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Useful walkers for traversing terms
**
**
**/

#include "utils/term_walkers.h"

#include "assert.h"
#include "utils/term_analysis.h"

using namespace smt;
using namespace std;

// include bv versions for boolector
// won't matter for other solvers because the sort won't be boolean
unordered_set<PrimOp> boolops(
    { And, Or, Xor, Not, Implies, BVAnd, BVOr, BVXor, BVNand, BVNot });

namespace pono {

void TermOpCollector::find_matching_terms(
    Term t, const unordered_set<PrimOp> & prim_ops, UnorderedTermSet & out)
{
  // set the pointers for use in visit/visit_term
  prim_ops_ = &prim_ops;
  out_ = &out;

  // visit all the subterms and collect them in out
  // see visit_term below
  visit(t);

  // set the pointers back to null
  // better to be null than stale
  prim_ops_ = nullptr;
  out_ = nullptr;
}

WalkerStepResult TermOpCollector::visit_term(smt::Term & term)
{
  // pointers should be non-null
  assert(prim_ops_);
  assert(out_);

  if (preorder_) {
    Op op = term->get_op();
    if (prim_ops_->find(op.prim_op) != prim_ops_->end()) {
      out_->insert(term);
    }
  }

  return Walker_Continue;
}

SubTermCollector::SubTermCollector(const smt::SmtSolver & solver,
                                   bool exclude_bools,
                                   bool exclude_funs,
                                   bool exclute_ites)
    : super(solver, true),
      exclude_bools_(exclude_bools),
      exclude_funs_(exclude_funs),
      exclude_ites_(exclute_ites),
      boolsort_(solver_->make_sort(BOOL))
{
}

void SubTermCollector::collect_subterms(Term term) { visit(term); }

WalkerStepResult SubTermCollector::visit_term(smt::Term & term)
{
  if (preorder_) {
    Sort sort = term->get_sort();

    // want to include literals as predicates here so passing true
    if (is_predicate(term, boolsort_, true)) {
      // TODO: consider handling ITEs here
      //       will still show up in predicates instead of being expanded
      predicates_.insert(term);
    } else if (exclude_bools_ && sort == boolsort_) {
      // special-case for boolector to make sure it keeps terms like
      // ((_ extract 0 0) x)
      // because it cannot distinguish between bit-vectors of width 1
      // and booleans

      bool add_term =
          !term->is_value();  // don't want to add constant true/false

      add_term &= (boolops.find(term->get_op().prim_op) == boolops.end());

      if (add_term) {
        // one extra case to consider is equalities between booleans
        // don't want to add those because it's just <->
        // e.g. it's technically a boolean operation
        Op op = term->get_op();
        if (op == Equal || op == BVComp) {
          Sort child_sort = (*(term->begin()))->get_sort();
          if (child_sort == boolsort_) {
            add_term = false;
          }
        }
      }

      if (add_term) {
        subterms_[sort].insert(term);
      }
    } else if ((exclude_funs_ && sort->get_sort_kind() == FUNCTION)
               || (exclude_ites_ && term->get_op() == Ite)) {
      return Walker_Continue;
    } else {
      subterms_[sort].insert(term);
    }
  }
  return Walker_Continue;
}

}  // namespace pono
