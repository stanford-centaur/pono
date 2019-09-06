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
    throw "Unknown symbols";
  }
  init_ = init;
  trans_ = trans;
}

void RelationalTransitionSystem::set_trans(const smt::Term trans)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(trans))
  {
    throw "Unknown symbols";
  }
  trans_ = trans;
}

void RelationalTransitionSystem::constrain_init(const smt::Term constraint)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(constraint))
  {
    throw "Unknown symbols";
  }
  init_ = solver_->make_term(And, init_, constraint);
}

void RelationalTransitionSystem::constrain_trans(const smt::Term constraint)
{
  // TODO: Only do this check in debug mode
  if (!known_symbols(constraint))
  {
    throw "Unknown symbols";
  }
  trans_ = solver_->make_term(And, trans_, constraint);
}

Term RelationalTransitionSystem::curr(const smt::Term term)
{
  UnorderedTermMap cache;
  TermVec to_visit {term};
  Term t;
  while(to_visit.size())
  {
    t = to_visit.back();
    to_visit.pop_back();
    if (cache.find(t) == cache.end())
    {
      // doesn't get updated yet, just marking as visited
      cache[t] = t;
      to_visit.push_back(t);
      for (auto c : t)
      {
        to_visit.push_back(c);
      }
    }
    else if(next_states_map_.find(t) != next_states_map_.end())
    {
      cache[t] = next_states_map_[t];
    }
    else
    {
      TermVec cached_children;
      for (auto c : t)
      {
        cached_children.push_back(cache.at(c));
      }
      Term rebuilt = solver_->make_term(t->get_op(), cached_children);
      cache[t] = rebuilt;
    }
  }

  return cache[term];
}

Term RelationalTransitionSystem::next(const smt::Term term)
{
  UnorderedTermMap cache;
  TermVec to_visit {term};
  Term t;
  while(to_visit.size())
  {
    t = to_visit.back();
    to_visit.pop_back();
    if (cache.find(t) == cache.end())
    {
      // doesn't get updated yet, just marking as visited
      cache[t] = t;
      to_visit.push_back(t);
      for (auto c : t)
      {
        to_visit.push_back(c);
      }
    }
    else if(states_map_.find(t) != states_map_.end())
    {
      cache[t] = states_map_[t];
    }
    else
    {
      TermVec cached_children;
      for (auto c : t)
      {
        cached_children.push_back(cache.at(c));
      }
      Term rebuilt = solver_->make_term(t->get_op(), cached_children);
      cache[t] = rebuilt;
    }
  }

  return cache[term];
}

bool RelationalTransitionSystem::is_curr_var(const smt::Term sv)
{
  return (states_.find(sv) != states_.end());
}

bool RelationalTransitionSystem::is_next_var(const smt::Term sv)
{
  return (next_states_.find(sv) != next_states_.end());
}

// overloaded methods (using next-state variables)
Term RelationalTransitionSystem::make_state(const string name, const Sort sort)
{
  Term state = solver_->make_term(name, sort);
  Term next_state = solver_->make_term(name + ".next", sort);
  states_.insert(state);
  next_states_.insert(next_state);
  states_map_[state] = next_state;
  next_states_map_[next_state] = state;
  return state;
}

bool RelationalTransitionSystem::known_symbols(const Term term)
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

// protected methods

bool RelationalTransitionSystem::no_next(smt::Term term)
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

}
