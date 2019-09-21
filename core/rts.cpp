#include "rts.h"

using namespace smt;
using namespace std;

namespace cosa
{

void RelationalTransitionSystem::set_behavior(const smt::Term init, const smt::Term trans)
{
  // TODO: Only do this check in debug mode
  if(!known_symbols(init) || !known_symbols(trans))
  {
    throw CosaException("Unknown symbols");
  }
  init_ = init;
  trans_ = trans;
}

void RelationalTransitionSystem::set_trans(const smt::Term trans)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(trans))
  {
    throw CosaException("Unknown symbols");
  }
  trans_ = trans;
}

void RelationalTransitionSystem::constrain_trans(const smt::Term constraint)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(constraint))
  {
    throw CosaException("Unknown symbols");
  }
  trans_ = solver_->make_term(And, trans_, constraint);
}

Term RelationalTransitionSystem::curr(const smt::Term term) const
{
  return solver_->substitute(term, curr_map_);
}

bool RelationalTransitionSystem::is_curr_var(const smt::Term sv) const
{
  return (states_.find(sv) != states_.end());
}

bool RelationalTransitionSystem::is_next_var(const smt::Term sv) const
{
  return (next_states_.find(sv) != next_states_.end());
}

// overloaded -- keep track of backwards mapping
Term RelationalTransitionSystem::make_input(const string name,
                                            const Sort sort) {
  Term input = solver_->make_term(name, sort);
  inputs_.insert(input);
  // for invariant constraints, need to assert over next inputs
  Term next_input = solver_->make_term(name + ".next", sort);
  next_map_[input] = next_input;
  curr_map_[next_input] = input;
  return input;
}

// overloaded methods (using next-state variables)
Term RelationalTransitionSystem::make_state(const string name, const Sort sort)
{
  Term state = solver_->make_term(name, sort);
  Term next_state = solver_->make_term(name + ".next", sort);
  states_.insert(state);
  next_states_.insert(next_state);
  next_map_[state] = next_state;
  curr_map_[next_state] = state;
  return state;
}

// protected methods

bool RelationalTransitionSystem::known_symbols(const Term term) const
{
  UnorderedTermSet visited;
  TermVec to_visit{term};
  Term t;
  while(to_visit.size())
  {
    t = to_visit.back();
    to_visit.pop_back();

    if(visited.find(t) != visited.end())
    {
      // cache hit
      continue;
    }

    if(t->is_symbolic_const() &&
       !((inputs_.find(t) != inputs_.end()) ||
         (states_.find(t) != states_.end()) ||
         (next_states_.find(t) != next_states_.end())
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

}
