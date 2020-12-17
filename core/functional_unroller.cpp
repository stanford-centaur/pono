/*********************                                                        */
/*! \file functional_unroller.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief An unroller implementation for functional transition systems that has
**        a configurable parameter for when to introduce new timed variables.
**
**
**/
#include "core/functional_unroller.h"

#include "assert.h"
#include "utils/exceptions.h"

using namespace smt;
using namespace std;

namespace pono {

FunctionalUnroller::FunctionalUnroller(const TransitionSystem & ts,
                                       const SmtSolver & solver,
                                       size_t interval)
    : super(ts, solver), interval_(interval)
{
  if (!ts.is_functional()) {
    throw PonoException(
        "Can only use FunctionalUnroller on a FunctionalTransitionSystem.");
  }
}

Term FunctionalUnroller::at_time(const Term & t, unsigned int k)
{
  if (!ts_.no_next(t)) {
    throw PonoException(
        "Functional unroller cannot replace next state variables");
  }
  return super::at_time(t, k);
}

UnorderedTermMap & FunctionalUnroller::var_cache_at_time(unsigned int k)
{
  const UnorderedTermMap & state_updates = ts_.state_updates();
  while (time_cache_.size() <= k) {
    time_cache_.push_back(UnorderedTermMap());
    UnorderedTermMap & subst = time_cache_.back();
    const unsigned int t = time_cache_.size() - 1;

    // create new state variables instead of substituting
    // every interval_ steps (if interval_ nonzero)
    bool create_new = (interval_ && (k % interval_ == 0));

    for (auto v : ts_.statevars()) {
      if (create_new) {
        Term new_v = var_at_time(v, t);
        subst[v] = new_v;
      } else {
        assert(t > 0);
        Term subst_v =
            solver_->substitute(state_updates.at(v), time_cache_.at(t - 1));
        subst[v] = subst_v;
      }
    }

    // always need to create new input variables
    for (auto v : ts_.inputvars()) {
      Term new_v = var_at_time(v, t);
      subst[v] = new_v;
    }
  }

  return time_cache_[k];
}

}  // namespace pono
