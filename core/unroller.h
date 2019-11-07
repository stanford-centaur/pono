#pragma once

#include "rts.h"

#include "smt-switch/smt.h"

namespace cosa {

/*
 * unroller inspired by the ic3ia unroller:
 * https://es-static.fbk.eu/people/griggio/ic3ia/
 *
 */

class Unroller
{
 public:
  Unroller(const RelationalTransitionSystem & ts, smt::SmtSolver & solver);
  ~Unroller();

  smt::Term at_time(const smt::Term & t, unsigned int k);

  smt::Term untime(const smt::Term & t) const;

 private:
  smt::Term var_at_time(const smt::Term & v, unsigned int k);
  smt::UnorderedTermMap & time_cache_at_time(unsigned int k);

  const RelationalTransitionSystem & ts_;
  smt::SmtSolver & solver_;

  typedef std::vector<smt::UnorderedTermMap> TimeCache;
  TimeCache time_cache_;
  TimeCache time_var_map_;
  smt::UnorderedTermMap untime_cache_;

};  // class Unroller

}  // namespace cosa
