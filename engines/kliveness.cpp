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

#include "kliveness.h"

#include <map>

#include "modifiers/mod_ts_prop.h"
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
  if (opt.promote_inputvars_) {
    ts_ = promote_inputvars(ts_);
    assert(!ts_.inputvars().empty());
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
    if (k < live_count_ + 1) {
      // k must be at least live_count_ + 1 to be able to
      // count live_count_ times
      logger.log(1,
                 "The given bound {} is too small to count {} times. "
                 "Please provide a larger bound.",
                 k,
                 live_count_);
      return ProverResult::UNKNOWN;
    }

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
    auto prop_and_monitor = instrument_ts(
        ts_k, tt.transfer_term(justice_conditions_.front()), live_count_);
    smt::Term safety_prop_term = prop_and_monitor.first;
    smt::Term monitor = prop_and_monitor.second;
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
    } else if (res == ProverResult::FALSE) {
      if (detect_revisit_in_cex(ts_k, safety_prover, monitor)) {
        return res;
      }
    }
  }
  return ProverResult::UNKNOWN;
}

std::pair<smt::Term, smt::Term> KLiveness::instrument_ts(
    TransitionSystem & ts_k,
    smt::Term justice_prop,
    const unsigned long & k) const
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
  return std::make_pair(safety_prop, counter);
}

std::shared_ptr<SafetyProver> KLiveness::make_safety_prover(
    const SafetyProperty & p, const TransitionSystem & ts) const
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

bool KLiveness::detect_revisit_in_cex(
    const TransitionSystem & ts,
    std::shared_ptr<SafetyProver> safety_prover,
    smt::Term counter) const
{
  if (live_count_ == 1) {
    // obviously no revisit when justice signal is triggered only once
    return false;
  }
  std::vector<smt::UnorderedTermMap> cex;
  safety_prover->witness(cex);
  assert(cex.size() >= live_count_ + 1);

  // identify the time steps where the justice signal is true
  std::vector<unsigned long> justice_steps;
  justice_steps.reserve(live_count_);
  smt::Term counter_val = cex.back().at(counter);
  for (std::size_t i = 1; i < cex.size(); ++i) {
    // traverse the cex backwards
    // if the counter value changes, it means the justice signal is true
    const std::size_t t = cex.size() - 1 - i;
    if (cex.at(t).at(counter) != counter_val) {
      justice_steps.push_back(t);
      counter_val = cex.at(t).at(counter);
    }
  }
  std::reverse(justice_steps.begin(), justice_steps.end());
  assert(justice_steps.size() == live_count_);

  // collect state values (excepts counter) at those time steps
  // and check if any two are identical
  std::map<smt::TermVec, std::size_t>
      visited;  // use map (to use unordered_map we have to
                // provide a hash function)
  smt::UnorderedTermSet state_vars;
  state_vars.reserve(ts.statevars().size()
                     - ts.statevars_with_no_update().size() - 1);
  for (const auto & sv : ts.statevars()) {
    if (sv != counter
        && ts.statevars_with_no_update().find(sv)
               == ts.statevars_with_no_update().end()) {
      state_vars.insert(sv);
    }
  }
  for (const auto & step : justice_steps) {
    smt::TermVec state_vals;
    state_vals.reserve(state_vars.size());
    for (const auto & sv : state_vars) {
      state_vals.push_back(cex.at(step).at(sv));
    }
    if (visited.find(state_vals) != visited.end()) {
      // found a revisit
      logger.log(2,
                 "Detected a lasso in the counterexample trace "
                 "(loop from step {} to {})",
                 visited.at(state_vals),
                 step);
      // TODO: store the lasso
#ifndef NDEBUG
      for (size_t t = 0; t <= step; ++t) {
        std::string step_info = (t == visited.at(state_vals)) ? " (loop start)"
                                : (t == step)                 ? " (revisit)"
                                                              : "";
        logger.log(3, "Step {}{}:", t, step_info);
        for (const auto & sv : state_vars) {
          logger.log(3, "  {} = {}", sv, cex.at(t).at(sv));
        }
      }
#endif
      return true;
    } else {
      visited.emplace(state_vals, step);
    }
  }
  return false;
}

}  // namespace pono
