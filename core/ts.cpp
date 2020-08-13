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

#include <functional>
#include "assert.h"

#include "ts.h"

using namespace smt;
using namespace std;

namespace pono {

TransitionSystem::TransitionSystem(const TransitionSystem & other_ts,
                                   TermTranslator & tt)
{
  function<Term(const Term &)> transfer;
  function<Term(const Term &, SortKind)> transfer_as;
  if (other_ts.solver() == tt.get_solver()) {
    // if the solvers are the same, don't need to transfer
    transfer = [](const Term & t) { return t; };
    // assume you don't need to do sort-casting for terms from the same solver
    transfer_as = [](const Term & t, SortKind sk) { return t; };
  } else {
    transfer = [&tt](const Term & t) { return tt.transfer_term(t); };
    transfer_as = [&tt](const Term & t, SortKind sk) {
      return tt.transfer_term(t, sk);
    };
  }

  solver_ = tt.get_solver();
  // transfer init and trans -- expect them to be boolean
  // will cast if underlying solver aliases Bool/BV1
  init_ = transfer_as(other_ts.init_, BOOL);
  trans_ = transfer_as(other_ts.trans_, BOOL);

  // populate data structures with translated terms

  for (auto v : other_ts.statevars_) {
    statevars_.insert(transfer(v));
  }

  for (auto v : other_ts.inputvars_) {
    inputvars_.insert(transfer(v));
  }

  for (auto v : other_ts.next_statevars_) {
    next_statevars_.insert(transfer(v));
  }

  for (auto elem : other_ts.named_terms_) {
    named_terms_[elem.first] = transfer(elem.second);
  }

  // variables might have already be in the TermTranslator cache
  // with a different sort (due to sort aliasing)
  // use the SortKind as a hint when transferring
  // sorts of the two terms should match for state updates and next_map
  Term key, val;
  for (auto elem : other_ts.state_updates_) {
    key = transfer(elem.first);
    val = transfer_as(elem.second, key->get_sort()->get_sort_kind());
    assert(key->get_sort() == val->get_sort());
    state_updates_[key] = val;
  }
  for (auto elem : other_ts.next_map_) {
    key = transfer(elem.first);
    val = transfer_as(elem.second, key->get_sort()->get_sort_kind());
    next_map_[key] = val;
  }

  for (auto elem : other_ts.curr_map_) {
    curr_map_[transfer(elem.first)] = transfer(elem.second);
  }

  functional_ = other_ts.functional_;
}

void TransitionSystem::set_init(const Term & init)
{
  // TODO: only do this check in debug mode
  if (!only_curr(init)) {
    throw PonoException(
        "Initial state constraints should only use current state variables");
  }

  init_ = init;
}

void TransitionSystem::constrain_init(const Term & constraint)
{
  // TODO: Only do this check in debug mode
  if (!only_curr(constraint)) {
    throw PonoException(
        "Initial state constraints should only use current state variables");
  }
  init_ = solver_->make_term(And, init_, constraint);
}

void TransitionSystem::assign_next(const Term & state, const Term & val)
{
  // TODO: only do this check in debug mode
  if (statevars_.find(state) == statevars_.end()) {
    throw PonoException("Unknown state variable");
  }

  if (!no_next(val)) {
    throw PonoException(
        "Got next state variable in RHS of functional assignment");
  }

  state_updates_[state] = val;
  trans_ = solver_->make_term(
      And, trans_, solver_->make_term(Equal, next_map_.at(state), val));
}

void TransitionSystem::add_invar(const Term & constraint)
{
  // TODO: only check this in debug mode
  if (only_curr(constraint)) {
    init_ = solver_->make_term(And, init_, constraint);
    trans_ = solver_->make_term(And, trans_, constraint);
    Term next_constraint = solver_->substitute(constraint, next_map_);
    // add the next-state version
    trans_ = solver_->make_term(And, trans_, next_constraint);
    constraints_.push_back(constraint);
    constraints_.push_back(next_constraint);
  } else {
    throw PonoException("Invariants should be over current states only.");
  }
}

void TransitionSystem::constrain_inputs(const Term & constraint)
{
  if (no_next(constraint)) {
    trans_ = solver_->make_term(And, trans_, constraint);
    constraints_.push_back(constraint);
  } else {
    throw PonoException("Cannot have next-states in an input constraint.");
  }
}

void TransitionSystem::add_constraint(const Term & constraint)
{
  if (only_curr(constraint)) {
    init_ = solver_->make_term(And, init_, constraint);
    trans_ = solver_->make_term(And, trans_, constraint);
    // add over next states
    Term next_constraint = solver_->substitute(constraint, next_map_);
    trans_ = solver_->make_term(And, trans_, next_constraint);
    constraints_.push_back(constraint);
    constraints_.push_back(next_constraint);
  } else if (no_next(constraint)) {
    trans_ = solver_->make_term(And, trans_, constraint);
    constraints_.push_back(constraint);
  } else {
    throw PonoException("Constraint cannot have next states");
  }
}

void TransitionSystem::name_term(const string name, const Term & t)
{
  if (named_terms_.find(name) != named_terms_.end()) {
    throw PonoException("Name " + name + " has already been used.");
  }
  named_terms_[name] = t;
}

Term TransitionSystem::make_inputvar(const string name, const Sort & sort)
{
  Term input = solver_->make_symbol(name, sort);
  inputvars_.insert(input);
  return input;
}

Term TransitionSystem::make_statevar(const string name, const Sort & sort)
{
  Term state = solver_->make_symbol(name, sort);
  Term next_state = solver_->make_symbol(name + ".next", sort);
  statevars_.insert(state);
  next_statevars_.insert(next_state);
  next_map_[state] = next_state;
  curr_map_[next_state] = state;
  return state;
}

Term TransitionSystem::curr(const Term & term) const
{
  return solver_->substitute(term, curr_map_);
}

Term TransitionSystem::next(const Term & term) const
{
  if (next_map_.find(term) != next_map_.end()) {
    return next_map_.at(term);
  }
  return solver_->substitute(term, next_map_);
}

bool TransitionSystem::is_curr_var(const Term & sv) const
{
  return (statevars_.find(sv) != statevars_.end());
}

bool TransitionSystem::is_next_var(const Term & sv) const
{
  return (next_statevars_.find(sv) != next_statevars_.end());
}

// term building methods -- forwards to SmtSolver solver_

Sort TransitionSystem::make_sort(const std::string name, uint64_t arity)
{
  return solver_->make_sort(name, arity);
}

Sort TransitionSystem::make_sort(const SortKind sk)
{
  return solver_->make_sort(sk);
}

Sort TransitionSystem::make_sort(const SortKind sk, uint64_t size)
{
  return solver_->make_sort(sk, size);
}

Sort TransitionSystem::make_sort(const SortKind sk, const Sort & sort1)
{
  return solver_->make_sort(sk, sort1);
}

Sort TransitionSystem::make_sort(const SortKind sk,
                                 const Sort & sort1,
                                 const Sort & sort2)
{
  return solver_->make_sort(sk, sort1, sort2);
}

Sort TransitionSystem::make_sort(const SortKind sk,
                                 const Sort & sort1,
                                 const Sort & sort2,
                                 const Sort & sort3)
{
  return solver_->make_sort(sk, sort1, sort2, sort3);
}

Sort TransitionSystem::make_sort(const SortKind sk, const SortVec & sorts)
{
  return solver_->make_sort(sk, sorts);
}

Term TransitionSystem::make_term(bool b) { return solver_->make_term(b); }

Term TransitionSystem::make_term(int64_t i, const Sort & sort)
{
  return solver_->make_term(i, sort);
}

Term TransitionSystem::make_term(const std::string val,
                                 const Sort & sort,
                                 uint64_t base)
{
  return solver_->make_term(val, sort, base);
}

Term TransitionSystem::make_term(const Term & val, const Sort & sort)
{
  return solver_->make_term(val, sort);
}

Term TransitionSystem::make_term(const Op op, const Term & t)
{
  return solver_->make_term(op, t);
}

Term TransitionSystem::make_term(const Op op, const Term & t0, const Term & t1)
{
  return solver_->make_term(op, t0, t1);
}

Term TransitionSystem::make_term(const Op op,
                                 const Term & t0,
                                 const Term & t1,
                                 const Term & t2)
{
  return solver_->make_term(op, t0, t1, t2);
}

Term TransitionSystem::make_term(const Op op, const TermVec & terms)
{
  return solver_->make_term(op, terms);
}

// protected methods

bool TransitionSystem::contains(const Term & term,
                                UnorderedTermSetPtrVec term_sets) const
{
  UnorderedTermSet visited;
  TermVec to_visit{ term };
  Term t;
  while (to_visit.size()) {
    t = to_visit.back();
    to_visit.pop_back();

    if (visited.find(t) != visited.end()) {
      // cache hit
      continue;
    }

    if (t->is_symbolic_const()) {
      bool in_atleast_one = false;
      for (auto ts : term_sets) {
        if (ts->find(t) != ts->end()) {
          in_atleast_one = true;
          break;
        }
      }

      if (!in_atleast_one) {
        return false;
      }
    }

    visited.insert(t);
    for (auto c : t) {
      to_visit.push_back(c);
    }
  }

  return true;
}

bool TransitionSystem::only_curr(const Term & term) const
{
  return contains(term, UnorderedTermSetPtrVec{ &statevars_ });
}

bool TransitionSystem::no_next(const Term & term) const
{
  return contains(term, UnorderedTermSetPtrVec{ &statevars_, &inputvars_ });
}

bool TransitionSystem::known_symbols(const Term & term) const
{
  return contains(
      term,
      UnorderedTermSetPtrVec{ &statevars_, &inputvars_, &next_statevars_ });
}

}  // namespace pono
