#include "unroller.h"

#include <assert.h>

using namespace smt;

namespace cosa {

Unroller::Unroller(const RelationalTransitionSystem & ts,
                   SmtSolver & solver)
    : ts_(ts), solver_(solver)
{
}

Unroller::~Unroller()
{
}

Term Unroller::at_time(const Term & t, unsigned int k)
{
  UnorderedTermMap & cache = time_cache_at_time(k);

  auto it = cache.find(t);
  if (it != cache.end())
  {
    return it->second;
  }

  Term ret = solver_->substitute(t, cache);
  untime_cache_[ret] = t;
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
  assert(time_var_map_.size() > k);

  UnorderedTermMap & cache = time_var_map_[k];
  auto it = cache.find(v);
  if (it != cache.end())
  {
    return it->second;
  }

  std::string name = v->to_string();
  name += "@" + std::to_string(k);
  Term timed_v = solver_->make_symbol(name, v->get_sort());
  cache[v] = timed_v;

  return timed_v;
}

UnorderedTermMap & Unroller::time_cache_at_time(unsigned int k)
{
  while (time_var_map_.size() <= k + 1)
  {
    time_var_map_.push_back(UnorderedTermMap());
  }

  while (time_cache_.size() <= k)
  {
    time_cache_.push_back(UnorderedTermMap());
    UnorderedTermMap & subst = time_cache_.back();

    if (subst.size() == 0)
    {
      for (auto v : ts_.states())
      {
        Term vn = ts_.next(v);
        Term new_v = var_at_time(v, k);
        Term new_vn = var_at_time(v, k + 1);
        subst[v] = new_v;
        subst[vn] = new_vn;
        untime_cache_[new_v] = v;
        untime_cache_[new_vn] = vn;
      }
      for (auto v : ts_.inputs())
      {
        Term vn = ts_.next(v);
        Term new_v = var_at_time(v, k);
        Term new_vn = var_at_time(v, k + 1);
        subst[v] = new_v;
        subst[vn] = new_vn;
        untime_cache_[new_v] = v;
        untime_cache_[new_vn] = vn;
      }
    }
  }

  return time_cache_[k];
}

}  // namespace cosa
