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

#include "ts.h"

using namespace smt;
using namespace std;

namespace cosa {

void TransitionSystem::set_init(const Term & init)
{
  // TODO: only do this check in debug mode
  if (!only_curr(init)) {
    throw CosaException(
        "Initial state constraints should only use current state variables");
  }

  init_ = init;
}

void TransitionSystem::constrain_init(const Term & constraint)
{
  // TODO: Only do this check in debug mode
  if (!only_curr(constraint)) {
    throw CosaException(
        "Initial state constraints should only use current state variables");
  }
  init_ = solver_->make_term(And, init_, constraint);
}

void TransitionSystem::assign_next(const Term & state, const Term & val)
{
  // TODO: only do this check in debug mode
  if (states_.find(state) == states_.end()) {
    throw CosaException("Unknown state variable");
  }

  if (!no_next(val)) {
    throw CosaException(
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
    // add the next-state version
    trans_ = solver_->make_term(
        And, trans_, solver_->substitute(constraint, next_map_));
  } else {
    throw CosaException("Invariants should be over current states only.");
  }
}

void TransitionSystem::constrain_inputs(const Term & constraint)
{
  if (no_next(constraint)) {
    trans_ = solver_->make_term(And, trans_, constraint);
  } else {
    throw CosaException("Cannot have next-states in an input constraint.");
  }
}

void TransitionSystem::add_constraint(const Term & constraint)
{
  if (only_curr(constraint)) {
    init_ = solver_->make_term(And, init_, constraint);
    trans_ = solver_->make_term(And, trans_, constraint);
    // add over next states
    trans_ = solver_->make_term(
        And, trans_, solver_->substitute(constraint, next_map_));
  } else if (no_next(constraint)) {
    trans_ = solver_->make_term(And, trans_, constraint);
  } else {
    throw CosaException("Constraint cannot have next states");
  }
}

void TransitionSystem::name_term(const string name, const Term & t)
{
  if (named_terms_.find(name) != named_terms_.end()) {
    throw CosaException("Name " + name + " has already been used.");
  }
  named_terms_[name] = t;
}

Term TransitionSystem::make_input(const string name, const Sort & sort)
{
  Term input = solver_->make_symbol(name, sort);
  inputs_.insert(input);
  return input;
}

Term TransitionSystem::make_state(const string name, const Sort & sort)
{
  Term state = solver_->make_symbol(name, sort);
  Term next_state = solver_->make_symbol(name + ".next", sort);
  states_.insert(state);
  next_states_.insert(next_state);
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
  return (states_.find(sv) != states_.end());
}

bool TransitionSystem::is_next_var(const Term & sv) const
{
  return (next_states_.find(sv) != next_states_.end());
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
  return contains(term, UnorderedTermSetPtrVec{ &states_ });
}

bool TransitionSystem::no_next(const Term & term) const
{
  return contains(term, UnorderedTermSetPtrVec{ &states_, &inputs_ });
}

bool TransitionSystem::known_symbols(const Term & term) const
{
  return contains(term,
                  UnorderedTermSetPtrVec{ &states_, &inputs_, &next_states_ });
}

}  // namespace cosa
