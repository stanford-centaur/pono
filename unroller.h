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

  smt::untime untime(smt::term t);

 private:

  const RelationalTransitionSystem &ts_;

  typedef std::vector<UnorderedTermMap *> TimeCache;
  TimeCache time_cache_;
  UnorderedTermMap untime_cache_;
  
}; // class Unroller

} // namespace cosa

#endif
