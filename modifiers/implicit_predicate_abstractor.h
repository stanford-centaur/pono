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

#include "smt-switch/utils.h"

#include "abstractor.h"

namespace pono {

class ImplicitPredicateAbstractor : public Abstractor
{
 public:
  ImplicitPredicateAbstractor(const TransitionSystem & conc_ts,
                              TransitionSystem & abs_ts,
                              smt::SmtSolver reducer_slv);

  smt::Term abstract(smt::Term & t) override;

  smt::Term concrete(smt::Term & t) override;

  /** Add a predicate to the abstraction
   *  @param pred the predicate to add (over concrete current state variables)
   *  @return the condition: pred(X') <-> pred(X^)
   *          that was added to the abstract transition relation
   *          this is for use in a procedure that wants to incrementally add
   *          predicates instead of re-adding the whole updated transition
   * relation
   */
  smt::Term add_predicate(const smt::Term & pred);

  /** Returns reference to vector of all current predicates over
   *  current state variables
   *  @return vector of predicates
   */
  const smt::TermVec & predicates() const { return predicates_; };

 protected:
  void do_abstraction() override;

  smt::Term predicate_refinement(const smt::Term & pred);

  const smt::SmtSolver & solver_;

  smt::UnsatCoreReducer reducer_;

  RelationalTransitionSystem & abs_rts_;
  ///< relational version of abs_ts_
  ///< this abstraction requires a relational system

  smt::TermVec predicates_;  ///< list of predicates in abstraction over current
                             ///< state vars
};

}  // namespace pono
