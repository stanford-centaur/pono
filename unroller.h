#ifndef UNROLLER_H
#define UNROLLER_H

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
  Unroller(RelationalTransitionSystem &ts);
  ~Unroller();

  smt::Term at_time(smt::Term t, unsigned int k);

  smt::Term untime(smt::Term t);

 private:

  smt::Term var_at_time(smt::Term v, unsigned int k);
  smt::UnorderedTermMap &time_cache_at_time(unsigned int k);

  RelationalTransitionSystem &ts_;
  smt::SmtSolver &solver_;
  
  typedef std::vector<smt::UnorderedTermMap *> TimeCache;
  TimeCache time_cache_;
  smt::UnorderedTermMap untime_cache_;
  
}; // class Unroller

} // namespace cosa

#endif
