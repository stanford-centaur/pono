/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#include <assert.h>
#include <algorithm>

#include "smt-switch/utils.h"

#include "unroller.h"

using namespace smt;
using namespace std;

namespace pono {

Unroller::Unroller(const TransitionSystem & ts, const string & time_identifier)
  : ts_(ts), solver_(ts.solver()), time_id_(time_identifier)
{
  num_vars_ = ts_.statevars().size();
  num_vars_ += ts_.inputvars().size();
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

  return solver_->substitute(t, cache);
}

Term Unroller::untime(const Term & t) const
{
  return solver_->substitute(t, untime_cache_);
}

size_t Unroller::get_var_time(const Term & v) const
{
  auto it = var_times_.find(v);
  if (it == var_times_.end()) {
    throw PonoException("Cannot get time of " + v->to_string());
  } else {
    return (*it).second;
  }
}

size_t Unroller::get_curr_time(const smt::Term & t) const
{
  UnorderedTermSet free_vars;
  get_free_symbolic_consts(t, free_vars);
  unordered_set<size_t> times;
  for (auto fv : free_vars) {
    assert(fv->is_symbolic_const());
    times.insert(get_var_time(fv));
  }
  size_t max = *std::max_element(times.begin(), times.end());
  size_t min = *std::min_element(times.begin(), times.end());
  if (max - min > 1) {
    throw PonoException("Cannot get current time of non-unrolled term.");
  }

  return min;
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
  name += time_id_ + std::to_string(k);
  Term timed_v = solver_->make_symbol(name, v->get_sort());
  cache[v] = timed_v;
  untime_cache_[timed_v] = v;
  var_times_[timed_v] = k;

  return timed_v;
}

UnorderedTermMap & Unroller::var_cache_at_time(unsigned int k)
{
  while (time_cache_.size() <= k) {
    time_cache_.push_back(UnorderedTermMap());
    UnorderedTermMap & subst = time_cache_.back();
    const unsigned int t = time_cache_.size() - 1;

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
  }

  UnorderedTermMap & var_cache = time_cache_[k];
  size_t current_num_vars = ts_.statevars().size();
  current_num_vars += ts_.inputvars().size();

  // if new variables are added to the transition system, we need to update the
  // timed-var-term-map
  if (current_num_vars > num_vars_) {
    num_vars_ = current_num_vars;
    size_t t = 0;
    for (UnorderedTermMap & st : time_cache_) {
      for (auto v : ts_.statevars()) {
        Term vn = ts_.next(v);
        Term new_v = var_at_time(v, t);
        Term new_vn = var_at_time(v, t + 1);
        st[v] = new_v;
        st[vn] = new_vn;
      }
      for (auto v : ts_.inputvars()) {
        Term new_v = var_at_time(v, t);
        st[v] = new_v;
      }

      ++t;
    }
  }

  return var_cache;
}

}  // namespace pono
