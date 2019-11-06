#pragma once

#include "prop.h"
#include "proverresult.h"
#include "rts.h"
#include "smt-switch/smt.h"
#include "unroller.h"

namespace cosa {
class Prover
{
 public:
  Prover(const Property & p, smt::SmtSolver & s);
  virtual ~Prover();

  virtual void initialize();

  virtual ProverResult check_until(int k) = 0;

  bool witness(std::vector<smt::UnorderedTermMap> & out);

  ProverResult prove();

 protected:
  const RelationalTransitionSystem & ts_;
  const Property & property_;

  smt::SmtSolver & solver_;
  Unroller unroller_;

  int reached_k_;

  smt::Term bad_;
};
}  // namespace cosa
