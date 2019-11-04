#include "fts.h"

using namespace smt;
using namespace std;

namespace cosa {

void FunctionalTransitionSystem::set_init(const Term & init)
{
  // TODO: only do this check in debug mode
  if (!known_symbols(init))
  {
    throw CosaException("Unknown symbols in formula");
  }

  init_ = init;
}

void FunctionalTransitionSystem::constrain_init(const Term & constraint)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(constraint))
  {
    throw CosaException("Unknown symbols");
  }
  init_ = solver_->make_term(And, init_, constraint);
}

void FunctionalTransitionSystem::set_next(const Term & state, const Term & val)
{
  // TODO: only do this check in debug mode
  if (states_.find(state) == states_.end())
  {
    throw CosaException("Unknown state variable");
  }

  // TODO: check only_curr on val in debug mode
  state_updates_[state] = val;
  trans_ = solver_->make_term(
      And, trans_, solver_->make_term(Equal, next_map_.at(state), val));
}

void FunctionalTransitionSystem::add_invar(const Term & constraint)
{
  // TODO: only check this in debug mode
  if (only_curr(constraint))
  {
    init_ = solver_->make_term(And, init_, constraint);
    trans_ = solver_->make_term(And, trans_, constraint);
    // add the next-state version
    trans_ = solver_->make_term(
        And, trans_, solver_->substitute(constraint, next_map_));
  }
  else
  {
    throw CosaException(
        "Invariants should be over current states and inputs only.");
  }
}

Term FunctionalTransitionSystem::make_input(const string name,
                                            const Sort & sort)
{
  Term input = solver_->make_symbol(name, sort);
  inputs_.insert(input);
  // for invariant constraints, need to assert over next inputs
  Term next_input = solver_->make_symbol(name + ".next", sort);
  next_map_[input] = next_input;
  return input;
}

Term FunctionalTransitionSystem::make_state(const string name,
                                            const Sort & sort)
{
  Term state = solver_->make_symbol(name, sort);
  Term next_state = solver_->make_symbol(name + ".next", sort);
  // this is never used, so it shouldn't hurt performance
  // only here for consistency with relational transition system states_ data
  // structure
  states_.insert(state);
  next_map_[state] = next_state;
  return state;
}

void FunctionalTransitionSystem::name_term(const string name, const Term & t)
{
  if (named_terms_.find(name) != named_terms_.end())
  {
    throw CosaException("Name has already been used.");
  }
  named_terms_[name] = t;
}

// protected methods

bool FunctionalTransitionSystem::only_curr(const Term & term) const
{
  UnorderedTermSet visited;
  TermVec to_visit{ term };
  Term t;
  while (to_visit.size())
  {
    t = to_visit.back();
    to_visit.pop_back();

    if (visited.find(t) != visited.end())
    {
      // cache hit
      continue;
    }

    if (t->is_symbolic_const() && (states_.find(t) == states_.end())
        && (inputs_.find(t) == inputs_.end()))
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

bool FunctionalTransitionSystem::known_symbols(const Term & term)
{
  UnorderedTermSet visited;
  TermVec to_visit{ term };
  Term t;
  while (to_visit.size())
  {
    t = to_visit.back();
    to_visit.pop_back();

    if (visited.find(term) != visited.end())
    {
      // cache hit
      continue;
    }

    if (t->is_symbolic_const()
        && !((inputs_.find(t) != inputs_.end())
             || (states_.find(t) != states_.end())))
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

Term FunctionalTransitionSystem::to_next_func(const Term & term)
{
  return solver_->substitute(term, state_updates_);
}

}  // namespace cosa
