#pragma once

#include <unordered_map>

#include "smt-switch/smt.h"

#include "fts.h"

namespace cosa
{

class RelationalTransitionSystem : public FunctionalTransitionSystem
{
 public:
   RelationalTransitionSystem(smt::SmtSolver & s)
     : FunctionalTransitionSystem(s) {}

  /* Sets init and trans to the provided values
   * @param init the new initial state constraints (boolean sort)
   * @param trans the new transition relation constraints (boolean sort)
  */
  void set_behavior(const smt::Term init, const smt::Term trans);

  /* Sets transition relation to the provided formula
   * @param trans the new transition relation
   */
  void set_trans(const smt::Term trans);

  /* Add to the transition relation constraints
   * @param constraint new constraint on transition relation
   */
  void constrain_trans(const smt::Term constraint);

  /* Map all next state variables to current state variables in the term
   * @param t the term to map
   * @return the term with all current state variables
   */
  smt::Term curr(const smt::Term term) const;

  /* @param sv the state variable to check
   * @return true if sv is a current state variable
   *
   * Returns false for any other term
   */
  bool is_curr_var(const smt::Term sv) const;

  /* @param sv the state variable to check
   * @return true if sv is a next state variable
   *
   * Returns false for any other term
   */
  bool is_next_var(const smt::Term sv) const;

  // getters

  /* Returns the transition relation
   * @return a boolean term representing the transition relation
   */
  smt::Term trans() const { return trans_; };

  // overloaded
  smt::Term make_input(const std::string name, const smt::Sort sort);

  // overloaded
  smt::Term make_state(const std::string name, const smt::Sort sort);

 protected:
  smt::UnorderedTermSet next_states_;
  // maps next back to curr
  smt::UnorderedTermMap curr_map_;

  // helpers and checkers

  // overloaded
  bool known_symbols(const smt::Term term) const;

};

}
