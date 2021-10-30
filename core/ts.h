/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
 ** This file is part of the pono project.
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

#include "smt-switch/cvc5_factory.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"

namespace pono {

class TransitionSystem
{
 public:
  /** use cvc5 by default (doesn't require logging so pass false)
   *  it supports the most theories and doesn't rewrite-on-the-fly or alias
   * sorts
   *  this makes it a great candidate for representing the TransitionSystem */
  TransitionSystem()
      : solver_(smt::Cvc5SolverFactory::create(false)),
        init_(solver_->make_term(true)),
        trans_(solver_->make_term(true)),
        functional_(false),
        deterministic_(false)
  {
  }

  TransitionSystem(const smt::SmtSolver & s)
      : solver_(s),
        init_(s->make_term(true)),
        trans_(s->make_term(true)),
        functional_(false),
        deterministic_(false)
  {
  }

  friend void swap(TransitionSystem & ts1, TransitionSystem & ts2);

  /** Copy assignment using
   *  copy-and-swap idiom
   */
  TransitionSystem & operator=(TransitionSystem other);

  /** Specialized copy-constructor that moves terms to new solver
   *  the solver is in the term translator
   *  the TransitinSystem doesn't keep the TermTranslator because
   *  it should not be connected to the other TransitionSystems
   *  i.e. it should not have any references to old Terms
   *  which would be kept in the TermTranslator cache
   *  @param other_ts the transition system to copy
   *  @param tt the term translator to use
   *  @return a transition system using the solver in tt
   */
  TransitionSystem(const TransitionSystem & other_ts, smt::TermTranslator & tt);

  virtual ~TransitionSystem(){};

  /** Equality comparison between TransitionSystems
   *  compares each member variable
   *  @param other the transition system to compare to
   *  @return true iff all member variables are equivalent
   */
  bool operator==(const TransitionSystem & other) const;

  /** Disquality comparison between TransitionSystems
   *  compares each member variable
   *  @param other the transition system to compare to
   *  @return true iff any member variables are different
   */
  bool operator!=(const TransitionSystem & other) const;

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
   * Throws a PonoException if:
   *  1) state is not a state variable
   *  2) val contains any next state variables (assign next is for functional
   * assignment)
   *  3) state has already been assigned a next state update
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
   *  @param constraint the constraint to add
   *  @param to_init_and_next whether it should be added to init and
   *         over next state variables as well
   *         (if it only contains state variables)
   * throws an exception if it has next states (should go in trans)
   */
  void add_constraint(const smt::Term & constraint,
                      bool to_init_and_next = true);

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
  smt::Term make_inputvar(const std::string name, const smt::Sort & sort);

  /* Create an state of a given sort
   * @param name the name of the state
   * @param sort the sort of the state
   * @return the current state variable
   *
   * Can get next state var with next(const smt::Term t)
   */
  smt::Term make_statevar(const std::string name, const smt::Sort & sort);

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

  /* @param t the term to check
   * @return true if t is an input variable
   *
   * Returns false for any other term
   */
  bool is_input_var(const smt::Term & t) const;

  /** Looks for a representative name for a term
   *  It searches for a name that was assigned to the term
   *  if it cannot be found, then it just returns the
   *  smt-lib to_string
   *  @param t the term to look for a name for
   *  @return a string for the term
   */
  std::string get_name(const smt::Term & t) const;

  /** Find a term by name in the transition system.
   *  searches current and next state variables, inputs,
   *  and named terms.
   *  Throws a PonoException if there is no matching term.
   *
   *  @param name the name to look for
   *  @return the matching term if found
   */
  smt::Term lookup(std::string name) const;

  /** Adds a state variable with the given next state
   *  @param cv the current state variable
   *  @param nv the next state variable
   */
  void add_statevar(const smt::Term & cv, const smt::Term & nv);

  /** Adds an input variable
   *  @param v the input variable
   */
  void add_inputvar(const smt::Term & v);

  // getters
  /* Returns const reference to solver */
  const smt::SmtSolver & solver() const { return solver_; };

  /* Gets a non-const reference to the solver */
  smt::SmtSolver & get_solver() { return solver_; };

  const smt::UnorderedTermSet & statevars() const { return statevars_; };

  const smt::UnorderedTermSet & inputvars() const { return inputvars_; };

  /* Returns the initial state constraints
   * @return a boolean term constraining the initial state
   */
  smt::Term init() const { return init_; };

  /* Returns the transition relation
   * @return a boolean term representing the transition relation
   */
  smt::Term trans() const { return trans_; };

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

  /** @return the constraints of the system
   *  Note: these do not include next-state variable updates or initial state
   * constraints
   *  Returned as a vector of pairs where for each element:
   *    first: is the constraint
   *    second: is a boolean saying whether it can be added over init / next
   * states This allows you to re-add the constraints by unpacking them and
   * passing to add_constraint, e.g. add_constraint(e.first, e.second)
   */
  const std::vector<std::pair<smt::Term, bool>> & constraints() const
  {
    return constraints_;
  };

  /** Whether the transition system is functional
   *  NOTE: This does *not* actually analyze the transition relation
   *  TODO possibly rename to be less confusing
   *  currently means state updates are always
   *    next_var := f(current_vars, inputs)
   *  however, it allows (certain) constraints still
   *  and does not require that every state has an update
   */
  bool is_functional() const { return functional_; };

  /** Whether the system is deterministic
   * this is a stronger condition than functional
   * TODO possibly rename to be less confusing
   *
   * deterministic (currently) means
   *   1) is_functional() is true
   *   2) every state has a next state for any input
   *      (i.e. no extra constraints)
   *   3) every state variable has an update function
   *       --> there exists exactly one next state
   *           if current vars and inputs are fixed
   */
  bool is_deterministic() const { return deterministic_; };

  /* Returns true iff all the symbols in the formula are current states */
  bool only_curr(const smt::Term & term) const;

  /* Returns true iff all the symbols in the formula are inputs and current
   * states */
  bool no_next(const smt::Term & term) const;

  /** EXPERTS ONLY
   *  Drop the state update for these variables and rebuild the system
   *  @param svs the state variables to drop updates for
   */
  void drop_state_updates(const smt::TermVec & svs);

  /** EXPERTS ONLY
   *  Turns an input variable into a state variable
   *  IMPORTANT: this does not retroactively change constraints
   *  e.g. if a constraint was not added to init because it
   *  contains an input variable
   *
   *  @param iv the input variable to promote
   *
   *  The input variable iv stays the same, but it will now
   *  be registered as a state variable.
   */
  void promote_inputvar(const smt::Term & iv);

  /** EXPERTS ONLY
   * Replace terms in the transition system with other terms
   *  Traverses all the data structures and updates them with
   *  the replacements
   *  Recommended usage is for experts only. Use this feature
   *    to cut out parts of the design by replacing terms with
   *    fresh inputs. Then cone of influence reduction can remove
   *    the unconnected parts of the transition system
   *  @param to_replace a mapping from terms in the transition
   *         system to their replacement.
   */
  void replace_terms(const smt::UnorderedTermMap & to_replace);

  // term building functionality -- forwards to the underlying SmtSolver
  // assumes all Terms/Sorts belong to solver_
  // e.g. should only use it on terms created by this TransitionSystem

  /* Make an uninterpreted sort
   * SMTLIB: (declare-sort <name> <arity>)
   * @param name the name of the sort
   * @param arity the arity of the sort
   * @return a Sort object
   */
  smt::Sort make_sort(const std::string name, uint64_t arity);

  /* Create a sort
   * @param sk the SortKind (BOOL, INT, REAL)
   * @return a Sort object
   */
  smt::Sort make_sort(const smt::SortKind sk);

  /* Create a sort
   * @param sk the SortKind (BV)
   * @param size (e.g. bitvector width for BV SortKind)
   * @return a Sort object
   */
  smt::Sort make_sort(const smt::SortKind sk, uint64_t size);

  /* Create a sort
   * @param sk the SortKind
   * @param sort1 first sort
   * @return a Sort object
   * this method is currently unused but kept for API consistency
   */
  smt::Sort make_sort(const smt::SortKind sk, const smt::Sort & sort1);

  /* Create a sort
   * @param sk the SortKind
   * @param sort1 first sort
   * @param sort2 second sort
   * @return a Sort object
   * When sk == ARRAY, sort1 is the index sort and sort2 is the element sort
   */
  smt::Sort make_sort(const smt::SortKind sk,
                      const smt::Sort & sort1,
                      const smt::Sort & sort2);

  /* Create a sort
   * @param sk the SortKind
   * @param sort1 first sort
   * @param sort2 second sort
   * @param sort3 third sort
   * @return a Sort object
   */
  smt::Sort make_sort(const smt::SortKind sk,
                      const smt::Sort & sort1,
                      const smt::Sort & sort2,
                      const smt::Sort & sort3);

  /* Create a sort
   * @param sk the SortKind (FUNCTION)
   * @param sorts a vector of sorts (for function SortKind, last sort is return
   * type)
   * @return a Sort object
   * Note: This is the only way to make a function sort
   */
  smt::Sort make_sort(const smt::SortKind sk, const smt::SortVec & sorts);

  /* Make a boolean value term
   * @param b boolean value
   * @return a value term with Sort BOOL and value b
   */
  smt::Term make_term(bool b);

  /* Make a bit-vector, int or real value term
   * @param i the value
   * @param sort the sort to create
   * @return a value term with Sort sort and value i
   */
  smt::Term make_term(int64_t i, const smt::Sort & sort);

  /* Make a bit-vector, int, real or (in the future) string value term
   * @param val the numeric value as a string, or a string value
   * @param sort the sort to create
   * @param base the base to interpret the value, for bit-vector sorts (ignored
   * otherwise)
   * @return a value term with Sort sort and value val
   */
  smt::Term make_term(const std::string val,
                      const smt::Sort & sort,
                      uint64_t base = 10);

  /* Make a value of a particular sort, such as constant arrays
   * @param val the Term used to create the value (.e.g constant array with 0)
   * @param sort the sort of value to create
   * @return a value term with Sort sort
   */
  smt::Term make_term(const smt::Term & val, const smt::Sort & sort);

  /* Make a new term
   * @param op the operator to use
   * @param t the child term
   * @return the created term
   */
  smt::Term make_term(const smt::Op op, const smt::Term & t);

  /* Make a new term
   * @param op the operator to use
   * @param t0 the first child term
   * @param t1 the second child term
   * @return the created term
   */
  smt::Term make_term(const smt::Op op,
                      const smt::Term & t0,
                      const smt::Term & t1);

  /* Make a new term
   * @param op the operator to use
   * @param t0 the first child term
   * @param t1 the second child term
   * @param t2 the third child term
   * @return the created term
   */
  smt::Term make_term(const smt::Op op,
                      const smt::Term & t0,
                      const smt::Term & t1,
                      const smt::Term & t2);

  /* Make a new term
   * @param op the operator to use
   * @param terms vector of children
   * @return the created term
   */
  smt::Term make_term(const smt::Op op, const smt::TermVec & terms);

  /* Rebuild transition relation 'trans_' based on set
     'state_vars_in_coi' of state-variables in cone-of-influence. The
     set 'state_vars_in_coi' is computed in the 'Prover' class that
     checks a property related to this transition system. Also, update
     the set of state/input variables to the passed sets
     'state_vars_in_coi' and 'input_vars_in_coi'. */
  void rebuild_trans_based_on_coi(
      const smt::UnorderedTermSet & state_vars_in_coi,
      const smt::UnorderedTermSet & input_vars_in_coi);

 protected:
  // solver
  smt::SmtSolver solver_;

  // initial state constraint
  smt::Term init_;

  // transition relation (functional in this class)
  smt::Term trans_;

  // system state variables
  smt::UnorderedTermSet statevars_;

  // set of next state variables
  smt::UnorderedTermSet next_statevars_;

  // system inputs
  smt::UnorderedTermSet inputvars_;

  // mapping from names to terms
  std::unordered_map<std::string, smt::Term> named_terms_;

  // mapping from terms to a representative name
  // because a term can have multiple names
  std::unordered_map<smt::Term, std::string> term_to_name_;

  // next state update function
  smt::UnorderedTermMap state_updates_;

  // maps states and inputs variables to next versions
  // note: the next state variables are only used
  //       on the left hand side of equalities in
  //       trans for functional transition systems
  smt::UnorderedTermMap next_map_;

  // maps next back to curr
  smt::UnorderedTermMap curr_map_;

  // whether the TransitionSystem is functional
  bool functional_;
  // whether the TransitionSystem is fully deterministic
  // i.e. if you fix the current states and inputs
  // is there only one next state
  // the only way functional_ could be true and deterministic_ be false
  // is if not all state variables have update functions
  bool deterministic_;

  // extra vector of terms to TransitionSystems that records constraints
  // added to the transition relation
  // For a functional system, you could now rebuild trans by AND-ing
  // together all the equalities from state_updates_
  // and these constraints
  // the boolean tells you whether it _can_ be added to next states
  // (the TransitionSystem will still check if it makes sense to add to
  // next states)
  // but this is crucial for some use cases: e.g. assuming the property
  // in the pre-state. It is very unsound to assume the property over
  // the init or next state variables, so the associated boolean for
  // that constraint would be false
  std::vector<std::pair<smt::Term, bool>> constraints_;
  ///< constraints added via
  ///< add_invar/constrain_inputs/add_constraint

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

  /* Returns true iff all the symbols in the formula are known */
  virtual bool known_symbols(const smt::Term & term) const;
};

}  // namespace pono
