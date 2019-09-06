#include "unroller.h"

namespace cosa {

  Unroller::Unroller(const RelationalTransitionSystem &ts):
    ts_(ts)
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
    smt::term ret;
    // TODO
    return ret;
  }

  smt::term Unroller::untime(smt::term t)
  {
    smt::term ret;
    // TODO
    return ret;
  }
  
} // namespace cosa
