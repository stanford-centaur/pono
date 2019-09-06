#include "fts.h"

using namespace smt;
using namespace std;

namespace cosa
{

void FunctionalTransitionSystem::set_init(const Term init)
{
  // TODO: only do this check in debug mode
  if (!known_symbols(init))
  {
    throw "Unknown symbols in formula";
  }

  init_ = init;
}

void FunctionalTransitionSystem::constrain_init(const smt::Term constraint)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(constraint))
  {
    throw "Unknown symbols";
  }
  init_ = solver_->make_term(And, init_, constraint);
}

void FunctionalTransitionSystem::set_next(const Term state, const Term val)
{
  // TODO: only do this check in debug mode
  if(states_.find(state) == states_.end())
  {
    throw "Unknown state variable";
  }

  state_updates_[state] = val;
}

void FunctionalTransitionSystem::add_constraint(const Term constraint)
{
  constraints_ = solver_->make_term(And, constraints_, constraint);
  // add the next-state version
  constraints_ = solver_->make_term(And, constraints_, to_next_func(constraint));
}

Term FunctionalTransitionSystem::make_input(const string name, const Sort sort)
{
  Term input = solver_->make_term(name, sort);
  inputs_.insert(input);
  return input;
}

Term FunctionalTransitionSystem::make_state(const string name, const Sort sort)
{
  Term state = solver_->make_term(name, sort);
  // this is never used, so it shouldn't hurt performance
  // only here for consistency with relational transition system states_ data structure
  states_.insert(state);
  return state;
}

void FunctionalTransitionSystem::name_term(const string name, const Term t)
{
  if (named_terms_.find(name) != named_terms_.end())
  {
    throw "Name has already been used.";
  }
  named_terms_[name] = t;
}

// protected methods

bool FunctionalTransitionSystem::known_symbols(const smt::Term term)
{
  UnorderedTermSet visited;
  TermVec to_visit{term};
  Term t;
  while(to_visit.size())
  {
    t = to_visit.back();
    to_visit.pop_back();

    if(visited.find(term) != visited.end())
    {
      // cache hit
      continue;
    }

    if(t->is_symbolic_const() &&
       !((inputs_.find(t) != inputs_.end()) ||
         (states_.find(t) != states_.end())
         ))
    {
      return false;
    }

    visited.insert(t);
    for (auto c : t)
    {
      to_visit.push_back(c);
    }
  }

  return true;
}

Term FunctionalTransitionSystem::to_next_func(Term term)
{
  return solver_->substitute(term, state_updates_);
}

}
