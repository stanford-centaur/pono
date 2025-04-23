#pragma once

#include "engines/prover.h"
#include "smt-switch/smt.h"

namespace pono {

class InterpSeqMC : public Prover
{
 public:
  InterpSeqMC(const Property & p,
              const TransitionSystem & ts,
              const smt::SmtSolver & slv,
              PonoOptions opt = PonoOptions());

  ~InterpSeqMC();

  typedef Prover super;

  void initialize() override;

  ProverResult check_until(int k) override;

 protected:
  bool step(int i);
  bool step_0();

  void reset_assertions(smt::SmtSolver & s);

  bool check_entail(const smt::Term & p, const smt::Term & q);

  smt::SmtSolver interpolator_;
  // for translating terms to interpolator_
  smt::TermTranslator to_interpolator_;
  // for translating terms to solver_
  smt::TermTranslator to_solver_;

  // set to true when a concrete_cex is found
  bool concrete_cex_;

  smt::Term init0_;
  smt::Term transA_;
  smt::Term transB_;
  smt::Term bad_disjuncts_;  ///< a disjunction of bads in the suffix

};  // class InterpSeqMC

}  // namespace pono
