#pragma once

#include "prop.h"
#include "prover.h"
#include "rts.h"
#include "unroller.h"

#include "proverresult.h"
#include "smt-switch/smt.h"

namespace cosa {

class InterpolationMC : public Prover
{
 public:
  // IMPORTANT: assume the property was built using the interpolating solver
  InterpolationMC(const Property & p,
                  smt::SmtSolver & itp,
                  smt::SmtSolver & slv);
  ~InterpolationMC();

  void initialize();

  ProverResult check_until(int k);

  ProverResult prove();

  bool witness(std::vector<smt::UnorderedTermMap> & out);

 private:
  const RelationalTransitionSystem & ts_;
  const Property & property_;

  bool step(int i);

  /* checks if the current Ri overapproximates R */
  bool check_overapprox();

  smt::SmtSolver & interpolator_;
  smt::SmtSolver & solver_;
  Unroller unroller_;

  int reached_k_;

  smt::UnorderedTermMap map_1_to_0;

  smt::Term init0_;
  smt::Term transA_;
  smt::Term transB_;
  smt::Term R_;
  smt::Term Ri_;

  smt::Term bad_;

};  // class InterpolationMC

}  // namespace cosa
