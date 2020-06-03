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
  UnorderedTermMap & cache = var_cache_at_time(k);

  // if t is a variable, it will be cached
  auto it = cache.find(t);
  if (it != cache.end()) {
    return it->second;
  }

  Term ret = solver_->substitute(t, cache);
  untime_cache_[ret] = t;

  return ret;
}

Term Unroller::untime(const Term & t) const
{
  return solver_->substitute(t, untime_cache_);
}

Term Unroller::var_at_time(const Term & v, unsigned int k)
{
  assert(v->is_symbolic_const());

  while (time_var_map_.size() <= k) {
    time_var_map_.push_back(UnorderedTermMap());
  }
  UnorderedTermMap & cache = time_var_map_[k];

  auto it = cache.find(v);
  if (it != cache.end()) {
    return it->second;
  }

  std::string name = v->to_string();
  name += "@" + std::to_string(k);
  Term timed_v = solver_->make_symbol(name, v->get_sort());
  cache[v] = timed_v;

  return timed_v;
}

UnorderedTermMap & Unroller::var_cache_at_time(unsigned int k)
{
  while (time_cache_.size() <= k) {
    time_cache_.push_back(UnorderedTermMap());
    UnorderedTermMap & subst = time_cache_.back();
    const unsigned int t = time_cache_.size() - 1;

    for (auto v : ts_.states()) {
      Term vn = ts_.next(v);
      Term new_v = var_at_time(v, t);
      Term new_vn = var_at_time(v, t + 1);
      subst[v] = new_v;
      subst[vn] = new_vn;
      untime_cache_[new_v] = v;
      untime_cache_[new_vn] = vn;
    }
    for (auto v : ts_.inputs()) {
      Term new_v = var_at_time(v, t);
      subst[v] = new_v;
      untime_cache_[new_v] = v;
    }
  }

  return time_cache_[k];
}

}  // namespace cosa
