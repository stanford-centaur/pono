#pragma once

#include "prop.h"
#include "prover.h"
#include "rts.h"
#include "unroller.h"

#include "proverresult.h"
#include "smt-switch/smt.h"

namespace cosa {

class InterpolantMC : public Prover
{
 public:
  // IMPORTANT: assume the property was built using the interpolating solver
  InterpolantMC(const Property & p, smt::SmtSolver & itp, smt::SmtSolver & slv);
  ~InterpolantMC();

  void initialize();

  ProverResult check_until(int k) override;
  bool witness(std::vector<smt::UnorderedTermMap> & out) override;

 private:
  const RelationalTransitionSystem & ts_;
  const Property & property_;

  bool step(int i);

  /* checks if the current Ri overapproximates R */
  bool check_overapprox();

  smt::SmtSolver & interpolator_;
  smt::SmtSolver & solver_;
  // for translating terms to solver_
  smt::TermTranslator tt_;
  Unroller unroller_;

  int reached_k_;

  // set to true when a concrete_cex is found
  bool concrete_cex_;

  smt::Term init0_;
  smt::Term transA_;
  smt::Term transB_;
  smt::Term R_;
  smt::Term Ri_;

  smt::Term bad_;

};  // class InterpolantMC

}  // namespace cosa
