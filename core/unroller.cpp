/*********************                                                        */
/*! \file 
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann
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


#include "unroller.h"

#include <assert.h>

using namespace smt;

namespace cosa {

Unroller::Unroller(const TransitionSystem & ts, SmtSolver & solver)
    : ts_(ts), solver_(solver)
{
}

Unroller::~Unroller() {}

Term Unroller::at_time(const Term & t, unsigned int k)
{
  UnorderedTermMap & cache = time_cache_at_time(k);

  auto it = cache.find(t);
  if (it != cache.end()) {
    return it->second;
  }

  Term ret = solver_->substitute(t, cache);
  cache[t] = ret;

  return ret;
}

Term Unroller::untime(const Term & t) const
{
  return solver_->substitute(t, untime_cache_);
}

Term Unroller::var_at_time(const Term & v, unsigned int k)
{
  assert(v->is_symbolic_const());

  UnorderedTermMap & cache = time_cache_[k];

  auto it = cache.find(v);
  if (it != cache.end()) {
    return it->second;
  }

  std::string name = v->to_string();
  name += "@" + std::to_string(k);
  Term timed_v = solver_->make_symbol(name, v->get_sort());
  cache[v] = timed_v;
  untime_cache_[timed_v] = v;

  return timed_v;
}

UnorderedTermMap & Unroller::time_cache_at_time(unsigned int k)
{
  size_t t = time_cache_.size();
  while (time_cache_.size() <= k+1) {
    time_cache_.push_back(UnorderedTermMap());
  }

  for (; t+1 < time_cache_.size(); ++t) {
    UnorderedTermMap & subst = time_cache_[t];

    for (auto &v : ts_.states()) {
      Term vn = ts_.next(v);
      Term new_v = var_at_time(v, t);
      Term new_vn = var_at_time(v, t + 1);
      subst[v] = new_v;
      subst[vn] = new_vn;
    }
    for (auto &v : ts_.inputs()) {
      Term new_v = var_at_time(v, t);
      subst[v] = new_v;
    }
  }

  return time_cache_[k];
}

}  // namespace cosa
