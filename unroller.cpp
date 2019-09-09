#include "unroller.h"

#include <assert.h>

namespace cosa {

  Unroller::Unroller(RelationalTransitionSystem &ts):
    ts_(ts),
    solver_(ts.solver())
  {
  }

  Unroller::~Unroller()
  {
    for (auto i : time_cache_) {
      delete i;
    }
  }

  smt::Term Unroller::at_time(smt::Term t, unsigned int k)
  {
    smt::UnorderedTermMap cache = time_cache_at_time(k);

    auto it = cache.find(t);
    if (it != cache.end()) {
      return it->second;
    }

    smt::Term ret = solver_->substitute(t, cache);
    untime_cache_[ret] = t;
    cache[t] = ret;

    return ret;
  }

  smt::Term Unroller::untime(smt::Term t)
  {
    return solver_->substitute(t, untime_cache_);
  }

  smt::Term Unroller::var_at_time(smt::Term v, unsigned int k)
  {
    assert(v->is_symbolic_const());

    std::string name = v->to_string();
    name += "@" + std::to_string(k);

    return solver_->make_term(name, v->get_sort());
  }

  smt::UnorderedTermMap &Unroller::time_cache_at_time(unsigned int k)
  {
    while (time_cache_.size() <= k) {
      time_cache_.push_back(new smt::UnorderedTermMap());
    }

    smt::UnorderedTermMap &subst = *(time_cache_[k]);

    if (subst.size() == 0) {
      const smt::UnorderedTermSet state_vars = ts_.states();
      const smt::UnorderedTermSet input_vars = ts_.inputs();
    
      for (auto it = state_vars.begin(), end = state_vars.end();
	   it != end; ++it) {
	smt::Term v = *it;
	smt::Term vn = ts_.next(v);
	smt::Term new_v = var_at_time(v, k);
	smt::Term new_vn = var_at_time(v, k+1);
	subst[v] = new_v;
	subst[vn] = new_vn;
	untime_cache_[new_v] = v;
	untime_cache_[new_vn] = vn;
      }
      for (auto it = input_vars.begin(), end = input_vars.end();
	   it != end; ++it) {
	smt::Term v = *it;
	smt::Term new_v = var_at_time(v, k);
	subst[v] = new_v;
	untime_cache_[new_v] = v;
      }
    }

    return subst;
  }
  
} // namespace cosa
