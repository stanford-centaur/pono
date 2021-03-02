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
                                       size_t interval,
                                       const string & time_identifier)
  : super(ts, time_identifier), interval_(interval),
    true_(solver_->make_term(true))
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
    extra_constraints_.push_back(true_);
    UnorderedTermMap & subst = time_cache_.back();
    const unsigned int t = time_cache_.size() - 1;
    assert(extra_constraints_.size() == t + 1);

    // create new state variables instead of substituting
    // every interval_ steps (if interval_ nonzero)
    // except when k is zero then we have to create_new regardless
    bool create_new = (interval_ && (k % interval_ == 0));
    create_new |= !t;

    for (auto v : ts_.statevars()) {
      bool no_update = state_updates.find(v) == state_updates.end();
      if (create_new || no_update) {
        Term new_v = var_at_time(v, t);
        subst[v] = new_v;
      }

      if (t == 0) {
        assert(create_new);  // should be creating new symbols at 0
        // no extra constraints at 0
        continue;
      } else if (no_update) {
        // nothing more to be done for implicit inputs
        continue;
      }

      assert(!no_update);
      Term fun_subst =
          solver_->substitute(state_updates.at(v), time_cache_.at(t - 1));

      if (create_new) {
        // add equality to extra constraints
        extra_constraints_[t] =
            solver_->make_term(And,
                               extra_constraints_[t],
                               solver_->make_term(Equal, subst[v], fun_subst));
      } else {
        assert(!subst[v]);  // expecting to not have been set already (e.g. be
                            // null)

        subst[v] = fun_subst;
      }
    }

    // always need to create new input variables
    for (auto v : ts_.inputvars()) {
      Term new_v = var_at_time(v, t);
      subst[v] = new_v;
    }
  }

  return time_cache_.at(k);
}

}  // namespace pono
