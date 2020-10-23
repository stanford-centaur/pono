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

namespace pono {

class ImplicitPredicateAbstractor : public Abstractor
{
 public:
  ImplicitPredicateAbstractor(const TransitionSystem & conc_ts,
                              TransitionSystem & abs_ts)
      : Abstractor(conc_ts, abs_ts), solver_(abs_ts.solver())
  {
    if (conc_ts_.solver() != abs_ts_.solver()) {
      throw PonoException(
          "For simplicity, expecting concrete and abstract system to use same "
          "solver.");
    }

    do_abstraction();
  }

  smt::Term abstract(smt::Term & t) override;

  smt::Term concrete(smt::Term & t) override;

  /** Add a predicate to the abstraction
   *  @param pred the predicate to add (over concrete current state variables)
   *  @return the condition: pred(X') <-> pred(X^)
   *          this is for use in a procedure that wants to incrementally add
   *          predicates instead of adding them all up front
   *  updates predicates_ and gets added as antecedent for predabs_label_
   */
  smt::Term add_predicate(const smt::Term & pred);

 protected:
  void do_abstraction() override;

  const smt::SmtSolver & solver_;

  smt::TermVec
      predicates_;  ///< list of predicates in abstraction over next state vars

  smt::Term
      predabs_label_;  ///< input variable in abs_ts_ that activates predicates
};

}  // namespace pono
