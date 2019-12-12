/*********************                                                        */
/*! \file 
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
 ** This file is part of the cosa2 project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief 
 **
 ** 
 **/


#pragma once

#include <string>
#include <unordered_map>

#include "smt-switch/smt.h"

#include "exceptions.h"

namespace cosa {

class FunctionalTransitionSystem
{
 public:
  FunctionalTransitionSystem(smt::SmtSolver & s)
      : solver_(s), init_(s->make_term(true)), trans_(s->make_term(true))
  {
  }

  /* Sets initial states to the provided formula
   * @param init the new initial state constraints
   */
  void set_init(const smt::Term & init);

  /* Add to the initial state constraints
   * @param constraint new constraint on initial states
   */
  void constrain_init(const smt::Term & constraint);

  /* Set the transition function of a state variable
   *   val is constrained to only use current state variables
   * @param state the state variable you are updating
   * @param val the value it should get
   */
  void set_next(const smt::Term & state, const smt::Term & val);

  /* Add constraint to the system
   * This is an invariant constraint, enforced over all time
   * @param constraint the boolean constraint term to add
   */
  void add_invar(const smt::Term & constraint);

  /* Create an input of a given sort
   * @param name the name of the input
   * @param sort the sort of the input
   * @return the input term
   */
  smt::Term make_input(const std::string name, const smt::Sort & sort);

  /* Create an state of a given sort
   * @param name the name of the state
   * @param sort the sort of the state
   * @return the current state variable
   *
   * Can get next state var with next(const smt::Term t)
   */
  smt::Term make_state(const std::string name, const smt::Sort & sort);

  /* Gives a term a name
   *   This can be used to track particular values in a witness
   * @param name the (unique) name to give the term
   * @param t the term to name
   *
   * Throws an exception if the name has already been used
   *  Note: giving multiple names to the same term is allowed
   */
  void name_term(const std::string name, const smt::Term & t);

  // getters
  smt::SmtSolver & solver() { return solver_; };

  const smt::UnorderedTermSet & states() const { return states_; };

  const smt::UnorderedTermSet & inputs() const { return inputs_; };

  /* Returns the initial state constraints
   * @return a boolean term constraining the initial state
   */
  smt::Term init() const { return init_; };

  /* Returns the next state updates
   * @return a map of functional next state updates
   */
  const smt::UnorderedTermMap & state_updates() const
  {
    return state_updates_;
  };

  /* @return the named terms mapping */
  const std::unordered_map<std::string, smt::Term> & named_terms() const
  {
    return named_terms_;
  };

 protected:
  // solver
  smt::SmtSolver & solver_;

  // initial state constraint
  smt::Term init_;

  // transition relation (functional in this class)
  smt::Term trans_;

  // next state update function
  smt::UnorderedTermMap state_updates_;

  // system state variables
  smt::UnorderedTermSet states_;

  // maps states and inputs variables to next versions
  // note: the next state variables are only used
  //       on the left hand side of equalities in
  //       trans for this functional class
  smt::UnorderedTermMap next_map_;

  // system inputs
  smt::UnorderedTermSet inputs_;

  // mapping from names to terms
  std::unordered_map<std::string, smt::Term> named_terms_;

  // helpers and checkers

  /* Returns true iff all the symbols in the formula are current states */
  bool only_curr(const smt::Term & term) const;

  /* Returns true iff all the symbols in the formula are known */
  bool known_symbols(const smt::Term & term);

  // Note: this method is not exposed, because the user should not be able to
  //       build terms with next states
  /* Replace all current states by their next-state updates, functionally */
  smt::Term to_next_func(const smt::Term & term);
};

}  // namespace cosa
