#pragma once

#include "rts.h"
#include "prop.h"
#include "smt-switch/smt.h"

#include "unroller.h"
#include "proverresult.h"

namespace cosa
{

class Bmc
{
public:
  Bmc(const Property &p, smt::SmtSolver &solver);
  ~Bmc();

  void initialize();
  
  ProverResult check_until(int k);

  ProverResult prove();
  
private:

  bool step(int i);
   
  const RelationalTransitionSystem &ts_;
  const Property &property_;

  smt::SmtSolver &solver_;

  Unroller unroller_;

  int reached_k_;

}; // class Bmc
  
} // namespace cosa


