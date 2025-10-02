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
 **
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
  smt::Term instrument_ts(TransitionSystem & ts_k,
                          smt::Term justice_prop,
                          const unsigned long & k);

  std::shared_ptr<SafetyProver> make_safety_prover(const SafetyProperty & p,
                                                   const TransitionSystem & ts);

  unsigned long live_count_;
};  // class KLiveness

}  // namespace pono
