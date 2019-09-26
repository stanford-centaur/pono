#pragma once

#include "prop.h"
#include "rts.h"
#include "unroller.h"

#include "proverresult.h"
#include "smt-switch/smt.h"

namespace cosa {
class Prover
{
 public:
  Prover(){};
  virtual ~Prover(){};

  virtual ProverResult check_until(int k) = 0;
  virtual ProverResult prove() = 0;
  virtual bool witness(std::vector<smt::UnorderedTermMap> & out) = 0;
};
}  // namespace cosa
