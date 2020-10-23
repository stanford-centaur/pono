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
      : Abstractor(conc_ts, abs_ts),
        solver_(abs_ts.solver()),
        abs_rts_(static_cast<RelationalTransitionSystem &>(abs_ts_))
  {
    if (conc_ts_.solver() != abs_ts_.solver()) {
      throw PonoException(
          "For simplicity, expecting concrete and abstract system to use same "
          "solver.");
    }

    // TODO: fix abstraction interface
    //       kind of strange to have to pass an empty abstract system

    if (abs_ts_.is_functional()) {
      throw PonoException(
          "Implicit predicate abstraction needs a relational abstract system");
    }

    // assume abstract transition starts empty
    // need to add all state variables and set behavior
    for (auto v : conc_ts_.statevars()) {
      abs_rts_.add_statevar(v, conc_ts_.next(v));
    }
    for (auto v : conc_ts_.inputvars()) {
      abs_rts_.add_inputvar(v);
    }
    // should start with the exact same behavior
    abs_rts_.set_behavior(conc_ts_.init(), conc_ts_.trans());

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

  /** Returns reference to vector of all current predicates over
   *  current state variables
   *  @return vector of predicates
   */
  const smt::TermVec & predicates() const { return predicates_; };

  /** Returns the input (more like definition) variable
   *  which activates the predicate equalities over next statevars
   *  and abstract next state vars
   *  e.g.
   *   predabs_label -> (p1(X) <-> p1(X^))
   *   predabs_label -> (p2(X) <-> p2(X^))
   *   etc...
   */
  const smt::Term & predabs_label() const { return predabs_label_; };

 protected:
  void do_abstraction() override;

  const smt::SmtSolver & solver_;

  RelationalTransitionSystem & abs_rts_;
  ///< relational version of abs_ts_
  ///< this abstraction requires a relational system

  smt::TermVec predicates_;  ///< list of predicates in abstraction over current
                             ///< state vars

  smt::Term predabs_label_;  ///< input variable in abs_ts_ that activates
                             ///< predicates over next state vars
};

}  // namespace pono
