/*********************                                                  */
/*! \file implicit_predicate_abstractor.h
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

#pragma once

#include "abstractor.h"
#include "core/unroller.h"
#include "smt-switch/term_translator.h"
#include "smt-switch/utils.h"

namespace pono {

class ImplicitPredicateAbstractor : public Abstractor
{
 public:
  ImplicitPredicateAbstractor(const TransitionSystem & conc_ts,
                              TransitionSystem & abs_ts,
                              Unroller & un);

  smt::Term abstract(smt::Term & t) override;

  smt::Term concrete(smt::Term & t) override;

  /** Returns the predicate refinement of the given predicate
   *  @param pred the predicate to refine
   *  (over concrete current state variables)
   *  @return the condition: pred(X') <-> pred(X^)
   *          this is for use in a procedure that wants to incrementally add
   *          predicates instead of re-adding the whole updated transition
   * relation
   */
  smt::Term predicate_refinement(const smt::Term & pred);

  void add_important_var(const smt::Term & v)
  {
    important_vars_.insert(v);
  }

  bool reduce_predicates(const smt::TermVec & cex,
                         const smt::TermVec & new_preds,
                         smt::TermVec & out);

  /** Does the abstraction and returns a set of concrete boolean symbols
   *  abstraction
   *  @return set of concrete boolean symbols
   */
  smt::UnorderedTermSet do_abstraction();

 protected:
  const smt::SmtSolver & solver_;

  Unroller & unroller_;

  smt::SmtSolver reducer_;
  smt::TermTranslator to_reducer_;

  RelationalTransitionSystem & abs_rts_;
  ///< relational version of abs_ts_
  ///< this abstraction requires a relational system

  bool abstracted_; ///< true iff do_abstraction has been called

  bool red_can_reset_;  ///< true iff reset_assertions workedo n reducer_

  smt::UnorderedTermSet important_vars_; ///< important variables
                                         ///< prioritize predicates containing these

  bool reset_reducer()
  {
    if (!red_can_reset_) {
      return false;
    }

    try {
      reducer_->reset_assertions();
      red_can_reset_ = true;
    }
    catch (SmtException & e) {
      red_can_reset_ = false;
    }
    return red_can_reset_;
  }
};

}  // namespace pono
