#ifndef FUNCTIONAL_TRANSITION_SYSTEM_H
#define FUNCTIONAL_TRANSITION_SYSTEM_H

#include <string>
#include <unordered_map>

#include "smt-switch/smt.h"

namespace cosa
{

class FunctionalTransitionSystem
{
 public:
    FunctionalTransitionSystem(smt::SmtSolver & s)
      : solver_(s),
        init_(s->make_value(true))
     {}

  /* Sets initial states to the provided formula
   * @param init the new initial state constraints
   */
  void set_init(const smt::Term init);

  /* Set the transition function of a state variable
   *   val is constrained to only use current state variables
   * @param state the state variable you are updating
   * @param val the value it should get
   */
  void set_next(const smt::Term state, const smt::Term val);

  /* Add constraint to the system
   * This is an invariant constraint, only over current states
   * @param constraint the boolean constraint term to add
  */
  void add_constraint(const smt::Term constraint);

  /* Create an input of a given sort
   * @param name the name of the input
   * @param sort the sort of the input
   * @return the input term
   */
  smt::Term make_input(const std::string name, const smt::Sort sort);

  /* Create an state of a given sort
   * @param name the name of the state
   * @param sort the sort of the state
   * @return the current state variable
   *
   * Can get next state var with next(const smt::Term t)
   */
  smt::Term make_state(const std::string name, const smt::Sort sort);

  /* Gives a term a name
   *   This can be used to track particular values in a witness
   * @param name the (unique) name to give the term
   * @param t the term to name
   *
   * Throws an exception if the name has already been used
   *  Note: giving multiple names to the same term is allowed
   */
  void name_term(const std::string name, const smt::Term t);

  // getters
  smt::SmtSolver & get_solver() { return solver_; };

  /* Returns the initial state constraints
   * @return a boolean term constraining the initial state
   */
  smt::Term init() { return init_; };

  /* Returns the next state updates
   * @return a map of functional next state updates
   */
  smt::UnorderedTermMap & state_updates() { return state_updates_; };

  /* @return the named terms mapping */
  std::unordered_map<std::string, smt::Term> & get_named_terms() { return named_terms_; };

 protected:
  smt::SmtSolver & solver_;
  smt::Term init_;
  smt::UnorderedTermMap state_updates_;
  smt::UnorderedTermSet states_;
  smt::UnorderedTermSet inputs_;
  std::unordered_map<std::string, smt::Term> named_terms_;

  // helpers and checkers

  /* Returns true iff all the symbols in the formula are known */
  bool known_symbols(const smt::Term term);

  // Note: this method is not exposed, because the user should not be able to
  //       build terms with next states
  /* Replace all current states by their next-state updates, functionally */
  smt::Term to_next_func(smt::Term term);

};

}

#endif
