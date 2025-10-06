/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Po-Chun Chien
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief Implementation of k-liveness
 **
 ** This liveness checking algorithm was introduced by Koen Claessen and Niklas
 ** Sörensson in FMCAD 2012 (https://ieeexplore.ieee.org/document/6462555).
 **/

#pragma once

#include "engines/prover.h"

namespace pono {

class KLiveness : public LivenessProver
{
 public:
  KLiveness(const LivenessProperty & p,
            const TransitionSystem & ts,
            const smt::SmtSolver & solver,
            PonoOptions opt = PonoOptions());

  ~KLiveness();

  typedef LivenessProver super;

  void initialize() override;

  ProverResult check_until(int k) override;

 protected:
  /** Instruments the original TS with counter.
   *  Returns a pair of <safety property, counter>.
   *  - Safety property: the counter never reaches the value k.
   *  - Counter: the introduced counter state variable.
   *  Note: the returned terms are in the context of the new TS (`ts_k`).
   */
  std::pair<smt::Term, smt::Term> instrument_ts(TransitionSystem & ts_k,
                                                smt::Term justice_prop,
                                                const unsigned long & k) const;

  std::shared_ptr<SafetyProver> make_safety_prover(
      const SafetyProperty & p, const TransitionSystem & ts) const;

  /** Detect if there is a lasso in the counterexample trace
   *  found by the safety prover
   */
  bool detect_revisit_in_cex(const TransitionSystem & ts,
                             std::shared_ptr<SafetyProver> safety_prover,
                             smt::Term counter) const;

  /* Perform BMC in lock-step to find a lasso */
  bool find_lasso_by_bmc(int bound) const;

  unsigned long live_count_;
  std::unique_ptr<SafetyProver> bmc_prover_;
};  // class KLiveness

}  // namespace pono
