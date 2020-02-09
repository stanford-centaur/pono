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


#include "fts.h"

using namespace smt;
using namespace std;

namespace cosa {

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

// protected methods

bool FunctionalTransitionSystem::known_symbols(const Term & term)
{
  UnorderedTermSet visited;
  TermVec to_visit{ term };
  Term t;
  while (to_visit.size()) {
    t = to_visit.back();
    to_visit.pop_back();

    if (visited.find(term) != visited.end()) {
      // cache hit
      continue;
    }

    if (t->is_symbolic_const()
        && !((inputs_.find(t) != inputs_.end())
             || (states_.find(t) != states_.end()))) {
      return false;
    }

    visited.insert(t);
    for (auto c : t) {
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
