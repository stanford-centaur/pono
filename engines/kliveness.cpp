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

#include "kliveness.h"

#include "modifiers/static_coi.h"
#include "utils/logger.h"
#include "utils/make_provers.h"

namespace pono {

KLiveness::KLiveness(const LivenessProperty & p,
                     const TransitionSystem & ts,
                     const smt::SmtSolver & solver,
                     PonoOptions opt)
    : super(p, ts, solver, opt)
{
  engine_ = opt.engine_;
  if (justice_conditions_.size() > 1) {
    throw PonoException("K-liveness only supports one justice condition.");
  }
  // copied from pono.cpp
  if (opt.static_coi_) {
    StaticConeOfInfluence coi(ts_, justice_conditions_, opt.verbosity_);
  }
}

KLiveness::~KLiveness() {}

void KLiveness::initialize()
{
  if (initialized_) {
    return;
  }
  super::initialize();
  live_count_ = 1;  // start with 1
}

ProverResult KLiveness::check_until(int k)
{
  initialize();
  for (; live_count_ <= options_.klive_bound_; ++live_count_) {
    logger.log(1, "k-liveness counter {}", live_count_);
    // To be extra safe, reset solver state here in case the safety prover does
    // not do it during initialization.
    solver_->reset_assertions();

    // instrument the TS with counter
    TransitionSystem ts_k = ts_;  // make a copy
    smt::Term safety_prop_term = instrument_ts(ts_k, live_count_);
    SafetyProperty safety_prop(solver_, safety_prop_term);

    // create a safety prover
    auto safety_prover = make_safety_prover(safety_prop, ts_k);
    assert(safety_prover);
    ProverResult res = safety_prover->check_until(k);
    // k-liveness can only confirm TRUE (unsat).
    // FALSE (sat) may be spurious and requires additional checks (not yet
    // implemented). For now, we simply continue with a higher live_count_.
    if (res == ProverResult::TRUE) {
      return res;
    }
  }
  return ProverResult::UNKNOWN;
}

smt::Term KLiveness::instrument_ts(TransitionSystem & ts_k,
                                   const unsigned long & k)
{
  throw PonoException("NYI");
}

std::shared_ptr<SafetyProver> KLiveness::make_safety_prover(
    const SafetyProperty & p, const TransitionSystem & ts)
{
  // modified from pono.cpp
  if (options_.cegp_abs_vals_) {
    return make_cegar_values_prover(engine_, p, ts, solver_, options_);
  } else if (options_.ceg_bv_arith_) {
    return make_cegar_bv_arith_prover(engine_, p, ts, solver_, options_);
  } else if (options_.ceg_prophecy_arrays_) {
    return make_ceg_proph_prover(engine_, p, ts, solver_, options_);
  } else {
    return make_prover(engine_, p, ts, solver_, options_);
  }
}

}  // namespace pono
