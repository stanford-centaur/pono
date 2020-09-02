/*********************                                                        */
/*! \file adaptive_unroller.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief An unroller implementation that will work even when new variables
**        are added to the transition system on the fly.
**
**/

#include "core/adaptive_unroller.h"

using namespace smt;
using namespace std;

namespace pono {

AdaptiveUnroller::AdaptiveUnroller(const TransitionSystem & ts,
                                   const smt::SmtSolver & solver)
    : super(ts, solver)
{
  num_vars_ = ts_.statevars().size();
  num_vars_ += ts_.inputvars().size();
}

UnorderedTermMap & AdaptiveUnroller::var_cache_at_time(unsigned int k)
{
  UnorderedTermMap & var_cache = super::var_cache_at_time(k);

  size_t current_num_vars = ts_.statevars().size();
  current_num_vars += ts_.inputvars().size();

  if (current_num_vars > num_vars_) {
    num_vars_ = current_num_vars;
    size_t t = 0;
    for (UnorderedTermMap & subst : time_cache_) {
      for (auto v : ts_.statevars()) {
        Term vn = ts_.next(v);
        Term new_v = var_at_time(v, t);
        Term new_vn = var_at_time(v, t + 1);
        subst[v] = new_v;
        subst[vn] = new_vn;
      }
      for (auto v : ts_.inputvars()) {
        Term new_v = var_at_time(v, t);
        subst[v] = new_v;
      }

      ++t;
    }
  }

  return var_cache;
}

}  // namespace pono
