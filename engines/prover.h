#pragma once

#include "proverresult.h"
#include "smt-switch/smt.h"

namespace cosa {
class Prover
{
 public:
  Prover();
  virtual ~Prover();

  virtual ProverResult check_until(int k) = 0;
  virtual bool witness(std::vector<smt::UnorderedTermMap> & out) = 0;

  ProverResult prove();
};
}  // namespace cosa
