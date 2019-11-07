#pragma once

#include "prover.h"

namespace cosa {

class Bmc : public Prover
{
 public:
  Bmc(const Property & p, smt::SmtSolver & solver);
  ~Bmc();

  typedef Prover super;

  void initialize() override;

  ProverResult check_until(int k) override;

 private:
  bool step(int i);

};  // class Bmc

}  // namespace cosa
