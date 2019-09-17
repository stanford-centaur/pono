#pragma once

#include "rts.h"
#include "prop.h"
#include "unroller.h"

#include "smt-switch/smt.h"
#include "proverresult.h"

namespace cosa
{
  class Prover
  {
  public:
    Prover() {};
    virtual ~Prover() {};

    virtual ProverResult check_until(int k) = 0;
    virtual ProverResult prove() = 0;
    virtual bool witness(std::vector<smt::UnorderedTermMap> &out) = 0;
  };
}
