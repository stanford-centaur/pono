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

using namespace smt;
using namespace std;

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

}  // namespace pono
