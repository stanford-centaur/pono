#include "unroller.h"

namespace cosa {

  Unroller::Unroller(const RelationalTransitionSystem &ts):
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

  smt::term Unroller::at_time(smt::term t, unsigned int k)
  {
    UnorderedTermMap &cache = time_cache_at_time(k);

    auto it = cache.find(t);
    if (it != cache.end()) {
      return it->second;
    }

    smt::term ret = solver_.substitute(t, *it);
    untime_cache_[ret] = t;
    cache[t] = ret;

    return ret;
  }

  smt::term Unroller::untime(smt::term t)
  {
    return solver_.substitute(t, untime_cache_);
  }

  smt::term Unroller::var_at_time(smt::term v, unsigned int k)
  {
    assert(t->is_symbolic_const());

    std::string name = t->to_string();
    name += "@" + std::to_string(k);

    return solver_.make_term(name, v->get_sort());
  }

  UnorderedTermMap &Unroller::time_cache_at_time(unsigned int k)
  {
    while (time_cache_.size() <= k) {
      time_cache_.push_back(new UnorderedTermMap());
    }

    UnorderedTermMap &subst = *(time_cache_[k]);

    if (subst.size() == 0) {
      const UnorderedTermSet state_vars = ts_.states();
      const UnorderedTermSet input_vars = ts_.inputs();
    
      for (size_t i = 0; i < state_vars.size(); ++i) {
	smt::term v = state_vars[i];
	smt::term vn = ts_.next(v);
	smt::term new_v = var_at_time(v, k);
	smt::term new_vn = var_at_time(v, k+1);
	subst[v] = new_v;
	subst[vn] = new_vn;
	untime_cache_[new_n] = v;
	untime_cache_[new_vn] = vn;
      }
      for (size_t i = 0; i < input_vars.size(); ++i) {
	smt::term v = input_vars[i];
	smt::term new_v = var_at_time(v, k);
	subst[v] = new_v;
	untime_cache_[new_v] = v;
      }
    }

    return subst;
  }
  
} // namespace cosa
