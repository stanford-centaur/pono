/*********************                                                        */
/*! \file term_walkers.h
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
#pragma once

#include "smt-switch/identity_walker.h"
#include "smt-switch/smt.h"

namespace pono {

/** Class for finding terms matching a given set of PrimOps
 *  uses protected inheritance because visit should not be
 *  called outside of this class (relies on the pointers
 *  being set)
 */
class TermOpCollector : protected smt::IdentityWalker
{
 public:
  TermOpCollector(const smt::SmtSolver & solver)
      // enable cache clearing because we want to traverse the whole
      // formula each time when looking for terms
      : smt::IdentityWalker(solver, true)
  {
  }

  /** Populates out with sub-terms of t that match one of the
   *  PrimOps in prim_ops
   *  @param t the term to traverse
   *  @param prim_ops the set of primitive ops to look for
   *  @param out the set to add matching subterms to
   */
  void find_matching_terms(smt::Term t,
                           const std::unordered_set<smt::PrimOp> & prim_ops,
                           smt::UnorderedTermSet & out);

 protected:
  smt::WalkerStepResult visit_term(smt::Term & term) override;

  const std::unordered_set<smt::PrimOp> * prim_ops_ = nullptr;
  smt::UnorderedTermSet * out_ = nullptr;
};

/** Class for collecting all subterms and grouping by sort
 *  It will also store predicates separately from all the other terms
 */
class SubTermCollector : public smt::IdentityWalker
{
 public:
  SubTermCollector(const smt::SmtSolver & solver,
                   bool exclude_bools = true,
                   bool exclude_funs = true,
                   bool exclude_ites = true);

  typedef smt::IdentityWalker super;

  void collect_subterms(smt::Term term);

  const std::unordered_map<smt::Sort, smt::UnorderedTermSet> & get_subterms() const
  {
    return subterms_;
  };

  const smt::UnorderedTermSet & get_predicates() const { return predicates_; }

 protected:
  bool exclude_bools_;  ///< if true, don't include boolean terms in subterms
                        ///<  (although predicates are still kept separately in
                        ///<  predicates_).

  bool exclude_funs_;  ///< if true, don't include function symbols

  bool exclude_ites_;  ///< if true, don't include ITEs

  smt::Sort boolsort_;  ///< boolean sort from solver_

  std::unordered_map<smt::Sort, smt::UnorderedTermSet> subterms_;

  smt::UnorderedTermSet predicates_;

  smt::WalkerStepResult visit_term(smt::Term & term) override;
};

}  // namespace pono
