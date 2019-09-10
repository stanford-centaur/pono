#pragma once

#include "rts.h"
#include "prop.h"
#include "smt-switch/smt.h"

#include "unroller.h"

namespace cosa
{

class Bmc
{
public:
  Bmc(const RelationalTransitionSystem &ts, const Property &p,
      smt::SmtSolver &solver);
  ~Bmc();

  bool check_until(size_t k);
  
private:

  bool step();
   
  const RelationalTransitionSystem &ts_;
  const Property &property_;

  smt::SmtSolver &solver_;

  Unroller unroller_;

}; // class Bmc
  
} // namespace cosa


