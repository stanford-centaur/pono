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

class TransitionSystem
{
 public:
  TransitionSystem(smt::SmtSolver & s)
      : solver_(s), init_(s->make_term(true)), trans_(s->make_term(true))
  {
  }

  virtual ~TransitionSystem(){};

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
   * Represents a functional update
   * @param state the state variable you are updating
   * @param val the value it should get
   */
  void assign_next(const smt::Term & state, const smt::Term & val);

  /* Add an invariant constraint to the system
   * This is enforced over all time
   * Specifically, it adds the constraint over both current and next variables
   * @param constraint the boolean constraint term to add (should contain only
   * state variables)
   */
  void add_invar(const smt::Term & constraint);

  /** Add a constraint over inputs
   * @param constraint to add (should not contain any next-state variables)
   */
  void constrain_inputs(const smt::Term & constraint);

  /** Convenience function for adding a constraint to the system
   * if the constraint only has current states, equivalent to add_invar
   * if there are only current states and inputs, equivalent to
   *   constrain_inputs
   * throws an exception if it has next states (should go in trans)
   */
  void add_constraint(const smt::Term & constraint);

  /* Gives a term a name
   *   This can be used to track particular values in a witness
   * @param name the (unique) name to give the term
   * @param t the term to name
   *
   * Throws an exception if the name has already been used
   *  Note: giving multiple names to the same term is allowed
   */
  void name_term(const std::string name, const smt::Term & t);

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

  /* Map all next state variables to current state variables in the term
   * @param t the term to map
   * @return the term with all current state variables
   */
  smt::Term curr(const smt::Term & term) const;

  /* Map all current state variables to next state variables in the term
   * @param t the term to map
   * @return the term with all next state variables
   */
  smt::Term next(const smt::Term & term) const;

  /* @param sv the state variable to check
   * @return true if sv is a current state variable
   *
   * Returns false for any other term
   */
  bool is_curr_var(const smt::Term & sv) const;

  /* @param sv the state variable to check
   * @return true if sv is a next state variable
   *
   * Returns false for any other term
   */
  bool is_next_var(const smt::Term & sv) const;

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

  /* Returns the transition relation
   * @return a boolean term representing the transition relation
   */
  smt::Term trans() const { return trans_; };

  /** Whether the transition system is functional
   *  NOTE: This does *not* actually analyze the transition relation
   *  it only returns true if it's a FunctionalTransitionSystem object
   */
  virtual bool is_functional() const { return false; };

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
  //       trans for functional transition systems
  smt::UnorderedTermMap next_map_;

  // system inputs
  smt::UnorderedTermSet inputs_;

  // mapping from names to terms
  std::unordered_map<std::string, smt::Term> named_terms_;
  smt::UnorderedTermSet next_states_;
  // maps next back to curr
  smt::UnorderedTermMap curr_map_;

  typedef std::vector<const smt::UnorderedTermSet *> UnorderedTermSetPtrVec;

  // helpers and checkers

  /** Returns true iff all symbols in term are present in at least one of the
   * term sets
   *  @param term the term to check
   *  @param term_sets a vector of sets to check membership for each symbol in
   * term
   *  @return true iff all symbols in term are in at least one of the term sets
   */
  bool contains(const smt::Term & term, UnorderedTermSetPtrVec term_sets) const;

  /* Returns true iff all the symbols in the formula are current states */
  bool only_curr(const smt::Term & term) const;

  /* Returns true iff all the symbols in the formula are inputs and current
   * states */
  bool no_next(const smt::Term & term) const;

  /* Returns true iff all the symbols in the formula are known */
  virtual bool known_symbols(const smt::Term & term) const;
};

}  // namespace cosa
