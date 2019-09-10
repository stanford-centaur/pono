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
  Bmc(const Property &p, smt::SmtSolver &solver);
  ~Bmc();

  void initialize();
  
  bool check_until(size_t k);

  bool prove();
  
private:

  bool step(size_t i);
   
  const RelationalTransitionSystem &ts_;
  const Property &property_;

  smt::SmtSolver &solver_;

  Unroller unroller_;

  size_t reached_k_;

}; // class Bmc
  
} // namespace cosa


