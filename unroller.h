#ifndef UNROLLER_H
#define UNROLLER_H

#include "rts.h"

#include "smt.h"

namespace cosa {
  
  /* 
   * unroller inspired by the ic3ia unroller:
   * https://es-static.fbk.eu/people/griggio/ic3ia/
   *
   */
  
class Unroller
{
 public:
  Unroller(const RelationalTransitionSystem &ts);

  smt::term at_time(smt::term t, unsigned int k);

  smt::term untime(smt::term t);

 private:

  smt::term var_at_time(smt::term v, unsigned int k);
  UnorderedTermMap time_cache_at_time(unsigned int k);

  const RelationalTransitionSystem &ts_;
  smt::SmtSolver &solver_;
  
  typedef std::vector<UnorderedTermMap *> TimeCache;
  TimeCache time_cache_;
  UnorderedTermMap untime_cache_;
  
}; // class Unroller

} // namespace cosa

#endif
