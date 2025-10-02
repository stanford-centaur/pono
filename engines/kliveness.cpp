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

#include <bit>

#include "modifiers/static_coi.h"
#include "smt/available_solvers.h"
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

    // create a new solver; provers may create named terms (e.g., unrolled state
    // vars) internally and could result in name conflicts
    smt::SmtSolver solver_k = create_solver_for(solver_->get_solver_enum(),
                                                engine_,
                                                false,
                                                options_.ceg_prophecy_arrays_,
                                                options_.printing_smt_solver_,
                                                options_.smt_solver_opts_);
    smt::TermTranslator tt(solver_k);

    // instrument the TS with counter
    TransitionSystem ts_k(ts_, tt);
    smt::Term safety_prop_term = instrument_ts(
        ts_k, tt.transfer_term(justice_conditions_.front()), live_count_);
    SafetyProperty safety_prop(solver_k, safety_prop_term);

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
                                   smt::Term justice_prop,
                                   const unsigned long & k)
{
  smt::Sort counter_sort;
  switch (options_.klive_counter_encoding_) {
    case KLivenessCounterEncoding::BV_BINARY:
      counter_sort = ts_k.make_sort(
          smt::SortKind::BV, (unsigned long)std::floor(std::log2(k)) + 1);
      break;
    case KLivenessCounterEncoding::BV_ONE_HOT:
      counter_sort = ts_k.make_sort(smt::SortKind::BV, k + 1);
      break;
    case KLivenessCounterEncoding::INTEGER:
      counter_sort = ts_k.make_sort(smt::SortKind::INT);
      break;
    default: throw PonoException("Unhandled k-liveness counter encoding");
  }

  smt::Term counter = ts_k.make_statevar(
      "pono_klive_counter_" + std::to_string(k), counter_sort);
  smt::Term zero = ts_k.make_term(0, counter_sort);
  smt::Term one = ts_k.make_term(1, counter_sort);

  smt::Term counter_init, counter_incr, count_to_val;
  switch (options_.klive_counter_encoding_) {
    case KLivenessCounterEncoding::BV_BINARY:
      counter_init = zero;
      counter_incr = ts_k.make_term(smt::PrimOp::BVAdd, counter, one);
      count_to_val = ts_k.make_term(k, counter_sort);
      break;
    case KLivenessCounterEncoding::BV_ONE_HOT:
      counter_init = one;
      counter_incr = ts_k.make_term(smt::PrimOp::BVShl, counter, one);
      count_to_val = ts_k.make_term(1UL << k, counter_sort);
      break;
    case KLivenessCounterEncoding::INTEGER:
      counter_init = zero;
      counter_incr = ts_k.make_term(smt::PrimOp::Plus, counter, one);
      count_to_val = ts_k.make_term(k, counter_sort);
      break;
    default: throw PonoException("Unhandled k-liveness counter encoding");
  }

  smt::Term counter_update =
      ts_k.make_term(smt::PrimOp::Ite, justice_prop, counter_incr, counter);
  ts_k.assign_next(counter, counter_update);
  ts_k.constrain_init(
      ts_k.make_term(smt::PrimOp::Equal, counter, counter_init));
  smt::Term safety_prop =
      ts_k.make_term(smt::PrimOp::Not,
                     ts_k.make_term(smt::PrimOp::Equal, counter, count_to_val));

  return safety_prop;
}

std::shared_ptr<SafetyProver> KLiveness::make_safety_prover(
    const SafetyProperty & p, const TransitionSystem & ts)
{
  // modified from pono.cpp
  if (options_.cegp_abs_vals_) {
    return make_cegar_values_prover(engine_, p, ts, ts.solver(), options_);
  } else if (options_.ceg_bv_arith_) {
    return make_cegar_bv_arith_prover(engine_, p, ts, ts.solver(), options_);
  } else if (options_.ceg_prophecy_arrays_) {
    return make_ceg_proph_prover(engine_, p, ts, ts.solver(), options_);
  } else {
    return make_prover(engine_, p, ts, ts.solver(), options_);
  }
}

}  // namespace pono
