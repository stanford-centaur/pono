#pragma once

#include "prop.h"
#include "rts.h"
#include "smt-switch/smt.h"

#include "prover.h"
#include "proverresult.h"
#include "unroller.h"

namespace cosa {

class Bmc : public Prover
{
 public:
  Bmc(const Property & p, smt::SmtSolver & solver);
  ~Bmc();

  void initialize();

  ProverResult check_until(int k) override;
  bool witness(std::vector<smt::UnorderedTermMap> & out) override;

 private:
  bool step(int i);

  const RelationalTransitionSystem & ts_;
  const Property & property_;

  smt::SmtSolver & solver_;

  Unroller unroller_;

  int reached_k_;

  smt::Term bad_;

};  // class Bmc

}  // namespace cosa
