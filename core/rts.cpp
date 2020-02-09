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


#include "rts.h"

using namespace smt;
using namespace std;

namespace cosa {

void RelationalTransitionSystem::set_behavior(const Term & init,
                                              const Term & trans)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(init) || !known_symbols(trans)) {
    throw CosaException("Unknown symbols");
  }
  init_ = init;
  trans_ = trans;
}

void RelationalTransitionSystem::set_init(const Term & init)
{
  // TODO: only do this check in debug mode
  if (!known_symbols(init)) {
    throw CosaException("Unknown symbols in formula");
  }

  init_ = init;
}

void RelationalTransitionSystem::constrain_init(const Term & constraint)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(constraint)) {
    throw CosaException("Unknown symbols");
  }
  init_ = solver_->make_term(And, init_, constraint);
}

void RelationalTransitionSystem::set_trans(const Term & trans)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(trans)) {
    throw CosaException("Unknown symbols");
  }
  trans_ = trans;
}

void RelationalTransitionSystem::constrain_trans(const Term & constraint)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(constraint)) {
    throw CosaException("Unknown symbols");
  }
  trans_ = solver_->make_term(And, trans_, constraint);
}

void RelationalTransitionSystem::set_next(const Term & state, const Term & val)
{
  // TODO: only do this check in debug mode
  if (states_.find(state) == states_.end()) {
    throw CosaException("Unknown state variable");
  }

  // TODO: check only_curr on val in debug mode
  state_updates_[state] = val;
  trans_ = solver_->make_term(
                              And, trans_, solver_->make_term(Equal, next_map_.at(state), val));
}

void RelationalTransitionSystem::add_invar(const Term & constraint)
{
  // TODO: only check this in debug mode
  if (only_curr(constraint)) {
    init_ = solver_->make_term(And, init_, constraint);
    trans_ = solver_->make_term(And, trans_, constraint);
    // add the next-state version
    trans_ = solver_->make_term(
                                And, trans_, solver_->substitute(constraint, next_map_));
  } else {
    throw CosaException(
                        "Invariants should be over current states and inputs only.");
  }
}

void RelationalTransitionSystem::name_term(const string name, const Term & t)
{
  if (named_terms_.find(name) != named_terms_.end()) {
    throw CosaException("Name has already been used.");
  }
  named_terms_[name] = t;
}

Term RelationalTransitionSystem::make_input(const string name,
                                            const Sort & sort)
{
  Term input = solver_->make_symbol(name, sort);
  inputs_.insert(input);
  // for invariant constraints, need to assert over next inputs
  Term next_input = solver_->make_symbol(name + ".next", sort);
  next_map_[input] = next_input;
  curr_map_[next_input] = input;
  return input;
}

Term RelationalTransitionSystem::make_state(const string name,
                                            const Sort & sort)
{
  Term state = solver_->make_symbol(name, sort);
  Term next_state = solver_->make_symbol(name + ".next", sort);
  states_.insert(state);
  next_states_.insert(next_state);
  next_map_[state] = next_state;
  curr_map_[next_state] = state;
  return state;
}

Term RelationalTransitionSystem::curr(const Term & term) const
{
  return solver_->substitute(term, curr_map_);
}

Term RelationalTransitionSystem::next(const Term & term) const
{
  if (next_map_.find(term) != next_map_.end()) {
    return next_map_.at(term);
  }
  return solver_->substitute(term, next_map_);
}

bool RelationalTransitionSystem::is_curr_var(const Term & sv) const
{
  return (states_.find(sv) != states_.end());
}

bool RelationalTransitionSystem::is_next_var(const Term & sv) const
{
  return (next_states_.find(sv) != next_states_.end());
}

// protected methods

bool RelationalTransitionSystem::only_curr(const Term & term) const
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

    if (t->is_symbolic_const() && (states_.find(t) == states_.end())
        && (inputs_.find(t) == inputs_.end())) {
      return false;
    }

    visited.insert(t);
    for (auto c : t) {
      to_visit.push_back(c);
    }
  }

  return true;
}

bool RelationalTransitionSystem::known_symbols(const Term & term) const
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

    if (t->is_symbolic_const()
        && !((inputs_.find(t) != inputs_.end())
             || (states_.find(t) != states_.end())
             || (next_states_.find(t) != next_states_.end()))) {
      return false;
    }

    visited.insert(t);
    for (auto c : t) {
      to_visit.push_back(c);
    }
  }

  return true;
}

}  // namespace cosa
