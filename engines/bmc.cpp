/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann, Florian Lonsing
 ** This file is part of the pono project.
 ** Copyright (c) 2019, 2021, 2022 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#include "bmc.h"

#include "utils/logger.h"

using namespace smt;

namespace pono {

Bmc::Bmc(const SafetyProperty & p,
         const TransitionSystem & ts,
         const SmtSolver & solver,
         PonoOptions opt)
    : super(p, ts, solver, opt)
{
  engine_ = Engine::BMC;
  bin_search_frames_ = 0;
  bound_step_ = opt.bmc_bound_step_;
  bound_start_ = opt.bmc_bound_start_;
}

Bmc::~Bmc() {}

void Bmc::initialize()
{
  if (initialized_) {
    return;
  }

  super::initialize();

  logger.log(2, "BMC adding init constraint for step 0");
  solver_->reset_assertions();
  solver_->assert_formula(unroller_.at_time(ts_.init(), 0));
}

void Bmc::reset_env()
{
  initialized_ = false;
  Bmc::initialize();
}

ProverResult Bmc::check_until(int k)
{
  initialize();

  // NOTE: There is a corner case where an instance is trivially unsatisfiable,
  // i.e., safe, when the conjunction of initial state predicate and transition
  // (+ any constraints) is already unsat. We could also check this using unsat
  // core functionality of solver (if supported), and check if bad state
  // predicate is in core.
  logger.log(1, "BMC bound_start_: {} ", bound_start_);
  logger.log(1, "BMC bound_step_: {} ", bound_step_);

  // Options 'bmc_exponential_step_' results in doubling the bound in every step
  const bool exp_step = options_.bmc_exponential_step_;

  for (int i = bound_start_; i <= k;
       i = exp_step ? (i == 0 ? 1 : i << 1) : (i + bound_step_)) {
    if (!step(i)) {
      compute_witness();
      return ProverResult::FALSE;
    }
  }
  return ProverResult::UNKNOWN;
}

bool Bmc::step(int i)
{
  logger.log(1, "\nBMC checking at bound: {}", i);

  if (i <= reached_k_) {
    return true;
  }

  bool res = true;
  if (i > 0) {
    // Add transitions depending on current interval '[reached_k_ + 1, i]'
    logger.log(2, "  BMC reached_k = {}, i = {} ", reached_k_, i);
    for (int j = reached_k_ == -1 ? 1 : reached_k_ + 1; j <= i; j++) {
      logger.log(2, "  BMC adding transition for j-1 = {}", j - 1);
      solver_->assert_formula(unroller_.at_time(ts_.trans(), j - 1));
      if (options_.bmc_neg_init_step_) {
        logger.log(2, "  BMC adding negated init constraint for step {}", j);
        Term not_init =
            solver_->make_term(PrimOp::Not, unroller_.at_time(ts_.init(), j));
        solver_->assert_formula(not_init);
      }
    }
  }

  solver_->push();

  // BMC is guaranteed to find a cex if we make sure to check all bounds
  // (default). This behavior can be overridden by
  // 'options_.bmc_single_bad_state_'
  const int cex_guarantee = !options_.bmc_single_bad_state_;

  Term clause;
  if (cex_guarantee) {
    // Make sure we cover *all* states by adding disjunctive bad state predicate
    clause = solver_->make_term(false);
    for (int j = reached_k_ + 1; j <= i; j++) {
      logger.log(2, "  BMC adding bad state constraint for j = {}", j);
      clause =
          solver_->make_term(PrimOp::Or, clause, unroller_.at_time(bad_, j));
    }
  } else {
    // Add a single bad state predicate (bugs might be missed)
    logger.log(2, "  BMC adding bad state constraint for i = {}", i);
    clause = unroller_.at_time(bad_, i);
  }

  solver_->assert_formula(clause);

  Result r = solver_->check_sat();
  if (r.is_sat()) {
    logger.log(1, "  BMC check at bound {} satisfiable", i);
    res = false;
    if (options_.bmc_allow_non_minimal_cex_) {
      // Terminate immediately; the reported bound 'reached_k' of the
      // cex is an upper bound of the actual cex within the interval
      // that was tested most recently. Option
      // 'options_.bmc_allow_non_minimal_cex' makes sense only if
      // 'options.bmc_bound_step > 1'.
      reached_k_ = i - 1;
      return res;
    }
    if (cex_guarantee) {
      logger.log(2, "  BMC saving reached_k_ = {}", reached_k_);
      int reached_k_saved = reached_k_;
      int cex_upper_bound = bmc_interval_get_cex_ub(reached_k_ + 1, i);
      // FIX for corner case of length-1 interval: some solvers don't
      // allow 'push' before 'get-value' even when no terms are asserted
      // before 'get-value'. Hence don't 'push' if we don't add any terms
      // in function 'bmc_interval_block_cex_ub'.
      if (cex_upper_bound + 1 <= i) solver_->push();
      bmc_interval_block_cex_ub(cex_upper_bound + 1, i);
      // Find shortest cex within tested interval given by bad state clause
      // 'success' will be false if binary search fails or linear search is
      // enabled
      bool success =
          options_.bmc_min_cex_linear_search_
              ? false
              : (!options_.bmc_min_cex_less_inc_bin_search_
                     ? find_shortest_cex_binary_search(cex_upper_bound)
                     : find_shortest_cex_binary_search_less_inc(
                           cex_upper_bound));
      if (!success) {
        assert(!options_.bmc_min_cex_less_inc_bin_search_);
        reached_k_ = reached_k_saved;
        // Clear constraints added during upper bound computation and binary
        // search
        if (cex_upper_bound + 1 <= i) solver_->pop();
        while (bin_search_frames_-- > 0) solver_->pop();
        find_shortest_cex_linear_search(cex_upper_bound);
      }
    } else {
      // Handle corner case when using single bad state constraints and interval
      // search: for witness printing, which depends on reached_k_, we must set
      // reached_k_ to the bound that precedes the bound 'i' where the cex was
      // found
      reached_k_ = i - 1;
    }
  } else {
    logger.log(1, "  BMC check at bound {} unsatisfiable", i);
    solver_->pop();
    // Optional: add negated bad state predicates to bounds where no cex was
    // found
    if (options_.bmc_neg_bad_step_ || options_.bmc_neg_bad_step_all_) {
      Term not_bad;
      if (options_.bmc_neg_bad_step_all_) {
        for (int j = reached_k_ + 1; j <= i; j++) {
          logger.log(
              2, "  BMC adding negated bad state constraint for j = {}", j);
          not_bad = solver_->make_term(PrimOp::Not, unroller_.at_time(bad_, j));
          solver_->assert_formula(not_bad);
        }
      } else {
        logger.log(
            2, "  BMC adding negated bad state constraint for i = {}", i);
        not_bad = solver_->make_term(PrimOp::Not, unroller_.at_time(bad_, i));
        solver_->assert_formula(not_bad);
      }
    }
    reached_k_ = i;
  }

  return res;
}

// Get an upper bound on the cex, which is located in interval '[lb,ub]'
int Bmc::bmc_interval_get_cex_ub(const int lb, const int ub)
{
  const Term true_term = solver_->make_term(true);
  assert(lb <= ub);

  logger.log(2,
             "  BMC get cex upper bound: lower bound = {}, upper bound = {} ",
             lb,
             ub);

  int j;
  for (j = lb; j <= ub; j++) {
    Term bad_state_at_j = unroller_.at_time(bad_, j);
    logger.log(2,
               "    BMC get cex upper bound, checking value of bad state "
               "constraint j = {}",
               j);
    if (solver_->get_value(bad_state_at_j) == true_term) {
      logger.log(2, "    BMC get cex upper bound, found at j = {}", j);
      break;
    }
  }
  assert(j <= ub);

  return j;
}

// Add negated bad state predicate for all bounds in interval '[start,end]'.
// This way, we restrict the search space of the solver to disregard these
// bounds when searching for a cex.
void Bmc::bmc_interval_block_cex_ub(const int start, const int end)
{
  logger.log(2,
             "  BMC permanently blocking interval [start,end] = [{},{}]",
             start,
             end);
  for (int k = start; k <= end; k++) {
    Term not_bad = solver_->make_term(PrimOp::Not, unroller_.at_time(bad_, k));
    logger.log(
        3,
        "    BMC adding permanent blocking bad state constraint for k = {}",
        k);
    solver_->assert_formula(not_bad);
  }
}

// Run binary search for cex within interval '[reached_k_ + 1, upper_bound]'.
bool Bmc::find_shortest_cex_binary_search(const int upper_bound)
{
  assert(bin_search_frames_ == 0);
  assert(reached_k_ < upper_bound);
  logger.log(2,
             "\n  BMC binary search, cex found in interval "
             "[reached_k+1,upper_bound] = [{},{}]",
             reached_k_ + 1,
             upper_bound);

  if (upper_bound - reached_k_ == 1) {
    logger.log(2,
               "  BMC interval has length 1, skipping search for shortest cex");
    return true;
  }

  int low = reached_k_ + 1;
  int high = upper_bound;
  while (low <= high) {
    int mid = low + (high - low) / 2;
    logger.log(2,
               "\n  BMC binary search, (low, mid, high) = ({}, {}, {})",
               low,
               mid,
               high);

    logger.log(3, "  BMC binary search, solver->push()");
    solver_->push();
    bin_search_frames_++;

    int j;
    // We search for cex in [low,mid] hence block [mid+1,high]
    logger.log(2,
               "  BMC binary search, searching for cex in [low,mid] = [{},{}]",
               low,
               mid);
    logger.log(
        2,
        "  BMC binary search, temporarily blocking [mid+1,high] = [{},{}]",
        mid + 1,
        high);
    for (j = mid + 1; j <= high; j++) {
      logger.log(3,
                 "  BMC binary search, finding shortest cex---"
                 "adding blocking bad state constraint for j = {}",
                 j);
      Term not_bad =
          solver_->make_term(PrimOp::Not, unroller_.at_time(bad_, j));
      solver_->assert_formula(not_bad);
    }

    Result r = solver_->check_sat();
    assert(r.is_sat() || r.is_unsat());
    if (r.is_sat()) {
      logger.log(2, "  BMC binary search, sat result: {}", r);
      logger.log(
          2, "  BMC binary search, cex found in [low,mid] = [{},{}]", low, mid);
      logger.log(
          2,
          "  BMC binary search, [mid+1,high] = [{},{}] now permanently blocked",
          mid + 1,
          high);
      // If low == mid in current iteration, then we have tested a single
      // bad state constraint and found a cex; can exit loop
      if (low == mid)
        break;
      else {
        high = bmc_interval_get_cex_ub(low, mid);
        bmc_interval_block_cex_ub(high + 1, mid);
      }
    } else {
      logger.log(2, "  BMC binary search, unsat result: {}", r);
      logger.log(
          2, "  BMC binary search, no cex in [low,mid] = [{},{}]", low, mid);
      if (low >= high) {
        // Handle rare corner case
        logger.log(1,
                   "  BMC binary search failure: formula overconstrained,"
                   " falling back to linear search");
        return false;
      }
      logger.log(2,
                 "  BMC binary search, unblocking [mid+1,high] = [{},{}]",
                 mid + 1,
                 high);
      // Remove previously added blocking constraints for [mid+1,high]
      logger.log(3, "  BMC binary search, solver->pop()");
      solver_->pop();
      assert(bin_search_frames_ > 0);
      bin_search_frames_--;
      // Update reached k; we have iteratively shown that no cex exists in
      // [0,mid]
      reached_k_ = mid;
      logger.log(
          2,
          "  BMC binary search, permanently blocking [low,mid] = [{},{}]",
          low,
          mid);
      // No cex found in [low,mid] hence block [low,mid]
      for (j = low; j <= mid; j++) {
        logger.log(3,
                   "  BMC binary search, finding shortest cex---"
                   "adding blocking bad state constraint for j = {}",
                   j);
        Term not_bad =
            solver_->make_term(PrimOp::Not, unroller_.at_time(bad_, j));
        solver_->assert_formula(not_bad);
      }

      low = mid + 1;
    }
  }

  // Must find cex inside sat-branch of if-then-else above
  assert(low <= high);
  // Reached_k_ has been correctly updated to low - 1, i.e., cex bound - 1
  assert(reached_k_ + 1 == low);
  logger.log(1,
             "  BMC binary search, shortest cex at bound low = {},"
             " reached_k = {}",
             low,
             reached_k_);
  return true;
}

/* Like 'find_shortest_cex_binary_search' but use new clause and new
   frame in each solver call to search for shortest
   counterexample. This has the effect that we potentially benefit less
   from incremental solving. */
bool Bmc::find_shortest_cex_binary_search_less_inc(const int upper_bound)
{
  assert(reached_k_ < upper_bound);
  logger.log(2,
             "  BMC less incremental binary search, cex found in interval "
             "[reached_k+1,upper_bound] = [{},{}]",
             reached_k_ + 1,
             upper_bound);

  if (upper_bound - reached_k_ == 1) {
    logger.log(2,
               "  BMC interval has length 1, skipping search for shortest cex");
    return true;
  }

  int low = reached_k_ + 1;
  int high = upper_bound;
  while (low <= high) {
    int mid = low + (high - low) / 2;
    logger.log(2,
               "\n  BMC binary search, (low, mid, high) = ({}, {}, {})",
               low,
               mid,
               high);

    logger.log(3, "  BMC binary search, solver->pop()");
    // Discard most recent bad state clause
    solver_->pop();
    logger.log(3, "  BMC binary search, solver->push()");
    solver_->push();

    int j;
    // We search for cex in [low,mid]
    logger.log(2,
               "  BMC binary search, searching for cex in [low,mid] = [{},{}]",
               low,
               mid);
    Term clause = solver_->make_term(false);
    for (j = low; j <= mid; j++) {
      logger.log(3,
                 "  BMC binary search, finding shortest cex---"
                 "adding bad state constraint for j = {}",
                 j);
      clause =
          solver_->make_term(PrimOp::Or, clause, unroller_.at_time(bad_, j));
    }
    solver_->assert_formula(clause);

    Result r = solver_->check_sat();
    assert(r.is_sat() || r.is_unsat());
    if (r.is_sat()) {
      logger.log(2, "  BMC binary search, sat result: {}", r);
      logger.log(
          2, "  BMC binary search, cex found in [low,mid] = [{},{}]", low, mid);
      // If low == mid in current iteration, then we have tested a single
      // bad state constraint; can exit loop in case of satisfiability
      if (low == mid)
        break;
      else {
        high = bmc_interval_get_cex_ub(low, mid);
      }
    } else {
      logger.log(2, "  BMC binary search, unsat result: {}", r);
      logger.log(
          2, "  BMC binary search, no cex in [low,mid] = [{},{}]", low, mid);
      // Update reached k; we have iteratively shown that no cex exists in
      // [0,mid]
      reached_k_ = mid;
      low = mid + 1;
    }
  }

  // Must find cex inside sat-branch of if-then-else above
  assert(low <= high);
  // 'reached_k_' has been correctly updated to low - 1, i.e., cex bound - 1
  assert(reached_k_ + 1 == low);
  logger.log(1,
             "  BMC binary search, shortest cex at bound low = {},"
             " reached_k = {}",
             low,
             reached_k_);
  return true;
}

// Run linear search for cex within interval '[reached_k_ + 1, upper_bound]'
void Bmc::find_shortest_cex_linear_search(const int upper_bound)
{
  assert(reached_k_ < upper_bound);
  logger.log(2,
             "  BMC linear search, cex found in interval: lower bound = "
             "reached k = {},"
             " upper bound = {}",
             reached_k_,
             upper_bound);

  if (upper_bound - reached_k_ == 1) {
    logger.log(2,
               "  BMC interval has length 1, skipping search for shortest cex");
    return;
  }

  int j;
  for (j = reached_k_ + 1; j <= upper_bound; j++) {
    // pop: remove the latest bad state clause
    solver_->pop();
    solver_->push();
    logger.log(
        2,
        "  BMC finding shortest cex---adding bad state constraint for j = {}",
        j);
    solver_->assert_formula(unroller_.at_time(bad_, j));
    if (solver_->check_sat().is_sat()) {
      break;
    } else
      reached_k_ = j;
  }
  // Must have found cex in the interval
  if (j > upper_bound)
    throw PonoException(
        "Unexpected BMC failure in linear search: formula overconstrained");

  assert(reached_k_ + 1 == j);
  logger.log(
      1,
      "  BMC linear search found shortest cex at bound j = {}, reached_k {}",
      j,
      reached_k_);
}

}  // namespace pono
