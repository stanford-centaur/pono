/*********************                                                  */
/*! \file interp_seq_mc.h
 ** \verbatim
 ** Top contributors (to current version):
 **   Po-Chun Chien
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** (in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief Implementation of Dual Approximated Reachability.
 **
 ** This unbounded model-checking algorithm was introduced by Yakir Vizel,
 ** Orna Grumberg, and Sharon Shoham in the paper "Intertwined Forward-Backward
 ** Reachability Analysis Using Interpolants" published in TACAS 2013
 ** (https://doi.org/10.1007/978-3-642-36742-7_22).
 **
 **/

#include "engines/dual_approx_reach.h"

#include "smt-switch/exceptions.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"

using namespace smt;

namespace pono {

DualApproxReach::DualApproxReach(const SafetyProperty & p,
                                 const TransitionSystem & ts,
                                 const SmtSolver & slv,
                                 PonoOptions opt)
    : super(p, ts, slv, opt),
      interpolator_(
          create_interpolating_solver_for(opt.smt_interpolator_, Engine::DAR)),
      to_interpolator_(interpolator_),
      to_solver_(solver_)
{
  engine_ = Engine::DAR;
}

DualApproxReach::~DualApproxReach() {}

void DualApproxReach::initialize()
{
  if (initialized_) {
    return;
  }

  super::initialize();

  interpolator_->reset_assertions();

  // register uninterpreted-function symbols,
  // (introduced during abstraction)
  UnorderedTermMap & cache = to_solver_.get_cache();
  UnorderedTermSet free_symbols;
  get_free_symbols(bad_, free_symbols);
  get_free_symbols(ts_.init(), free_symbols);
  get_free_symbols(ts_.trans(), free_symbols);
  for (const auto & s : free_symbols) {
    if (s->get_sort()->get_sort_kind() == FUNCTION) {
      cache[to_interpolator_.transfer_term(s)] = s;
    }
  }

  concrete_cex_ = false;
  forward_seq_.push_back(ts_.init());
  backward_seq_.push_back(bad_);
}

ProverResult DualApproxReach::check_until(int k)
{
  initialize();

  try {
    for (int i = 0; i <= k; ++i) {
      if (step(i)) {
        return ProverResult::TRUE;
      } else if (concrete_cex_) {
        compute_witness();
        return ProverResult::FALSE;
      }
    }
  }
  catch (InternalSolverException & e) {
    logger.log(1, "DAR failed due to: {}", e.what());
  }
  return ProverResult::UNKNOWN;
}

bool DualApproxReach::step(int i)
{
  if (i <= reached_k_) {
    return false;
  }
  assert(i == reached_k_ + 1);

  logger.log(1, "Running DAR at bound: {}", i);

  if (i == 0) {
    // Can't get an interpolant at bound 0
    // only checking for trivial bug
    return step_0();
  }

  update_term_map(i);
  if (!local_strengthen() && (i == 1 || !global_strengthen())) {
    concrete_cex_ = true;
    return false;
  }

  ++reached_k_;
  assert(i == reached_k_);
  assert(forward_seq_.size() == reached_k_ + 1);
  assert(forward_seq_.size() == backward_seq_.size());

  return check_fixed_point(forward_seq_) || check_fixed_point(backward_seq_);
}

bool DualApproxReach::step_0()
{
  solver_->reset_assertions();
  // push the unrolled formulas here
  // as compute_witness() rely on timed variables
  solver_->assert_formula(unroller_.at_time(forward_seq_.at(0), 0));   // init
  solver_->assert_formula(unroller_.at_time(backward_seq_.at(0), 0));  // bad

  Result r = solver_->check_sat();
  if (r.is_unsat()) {
    reached_k_ = 0;
  } else {
    concrete_cex_ = true;
  }
  return false;
}

void DualApproxReach::update_term_map(size_t i)
{
  // symbols are already created in solver
  // need to add symbols at the given time step to cache
  UnorderedTermMap & cache = to_solver_.get_cache();
  Term term;
  for (const auto & sv : ts_.statevars()) {
    term = unroller_.at_time(sv, i);
    cache[to_interpolator_.transfer_term(term)] = term;
  }
  for (const auto & iv : ts_.inputvars()) {
    term = unroller_.at_time(iv, i);
    cache[to_interpolator_.transfer_term(term)] = term;
  }
}

bool DualApproxReach::local_strengthen()
{
  // We want to find an index such that
  // forward_seq_[len-1-i](s0) & TR(s0, s1) & backward_seq_[i](s1) is unsat.
  // The search can be done in arbitrary order.
  // Here we follow the order from the original paper,
  // starting from the last element in forward_seq_
  size_t unsat_idx = 0;
  const size_t seq_len = forward_seq_.size();
  for (; unsat_idx < seq_len; ++unsat_idx) {
    solver_->reset_assertions();
    Term f0 = unroller_.at_time(forward_seq_.at(seq_len - 1 - unsat_idx), 0);
    Term b1 = unroller_.at_time(backward_seq_.at(unsat_idx), 1);
    Term tr = unroller_.at_time(ts_.trans(), 0);
    solver_->assert_formula(solver_->make_term(And, f0, tr, b1));
    Result r = solver_->check_sat();
    if (r.is_unsat()) {
      // found an index such that the conjunction is unsat
      break;
    }
  }
  if (unsat_idx == seq_len) {
    // no such index found, return false
    logger.log(1, "DAR: local strengthening failed");
    return false;
  }
  pairwise_strengthen(seq_len - 1 - unsat_idx);
  return true;
}

bool DualApproxReach::global_strengthen()
{
  assert(forward_seq_.size());
  TermVec int_assertions;
  int_assertions.reserve(forward_seq_.size() + 1);
  solver_->reset_assertions();
  solver_->push();
  Term unrolled_trans = unroller_.at_time(
      solver_->make_term(And, forward_seq_.at(0), ts_.trans()), 0);
  int_assertions.push_back(to_interpolator_.transfer_term(unrolled_trans));
  solver_->assert_formula(unrolled_trans);  // init(0) & TR(0, 1)

  const size_t seq_len = forward_seq_.size();
  size_t unsat_idx = 1;
  for (; unsat_idx < seq_len; ++unsat_idx) {
    solver_->push();
    unrolled_trans = unroller_.at_time(ts_.trans(), unsat_idx);
    solver_->assert_formula(unrolled_trans);
    int_assertions.push_back(to_interpolator_.transfer_term(unrolled_trans));
    solver_->push();
    Term ub = unroller_.at_time(backward_seq_.at(seq_len - 1 - unsat_idx),
                                unsat_idx + 1);
    solver_->assert_formula(ub);
    Result r = solver_->check_sat();
    if (r.is_unsat()) {
      // found an index for refinement
      int_assertions.push_back(to_interpolator_.transfer_term(ub));
      break;
    }
    if (unsat_idx == seq_len - 1) {
      // this is a concrete counterexample
      logger.log(1, "DAR: Found a concrete CEX");
      return false;
    }
    solver_->pop();  // pop backward_seq_ assertion
  }
  TermVec int_itpseq;
  Result int_r =
      interpolator_->get_sequence_interpolants(int_assertions, int_itpseq);
  if (!int_r.is_unsat()) {
    throw PonoException(
        "DAR: global strengthening failed, expect UNSAT interpolation query");
  }
  for (size_t i = 1; i < std::min(seq_len, unsat_idx + 1); ++i) {
    Term int_itp = int_itpseq.at(i - 1);
    Term itp = to_solver_.transfer_term(int_itp);
    logger.log(3,
               "DAR: strengthening forward reachability sequence at position "
               "{} with {}",
               i,
               itp);
    forward_seq_.at(i) =
        solver_->make_term(And, forward_seq_.at(i), unroller_.untime(itp));
  }
  pairwise_strengthen(unsat_idx);
  return true;
}

void DualApproxReach::pairwise_strengthen(const size_t idx)
{
  assert(forward_seq_.size() == backward_seq_.size());
  assert(idx < forward_seq_.size());
  const size_t seq_len = forward_seq_.size();
  Term tr = unroller_.at_time(ts_.trans(), 0);
  Term int_tr = to_interpolator_.transfer_term(tr);

  // strengthen and extend forward_seq_
  forward_seq_.push_back(solver_->make_term(true));
  for (size_t i = idx; i < seq_len; ++i) {
    Term f = unroller_.at_time(forward_seq_.at(i), 0);
    Term int_f = to_interpolator_.transfer_term(f);
    Term b = unroller_.at_time(backward_seq_.at(seq_len - 1 - i), 1);
    Term int_b = to_interpolator_.transfer_term(b);
    Term int_itp;
    Result r = interpolator_->get_interpolant(
        interpolator_->make_term(And, int_f, int_tr), int_b, int_itp);
    if (!r.is_unsat()) {
      throw PonoException(
          "DAR: pairwise strengthening failed, "
          "expect UNSAT interpolation query");
    }
    Term itp = to_solver_.transfer_term(int_itp);
    logger.log(3,
               "DAR: strengthening forward reachability sequence at position "
               "{} with {}",
               i + 1,
               itp);
    forward_seq_.at(i + 1) =
        solver_->make_term(And, forward_seq_.at(i + 1), unroller_.untime(itp));
  }

  // strengthen and extend backward_seq_
  backward_seq_.push_back(solver_->make_term(true));
  for (size_t i = seq_len - 1 - idx; i < seq_len; ++i) {
    Term f = unroller_.at_time(forward_seq_.at(seq_len - 1 - i), 0);
    Term int_f = to_interpolator_.transfer_term(f);
    Term b = unroller_.at_time(backward_seq_.at(i), 1);
    Term int_b = to_interpolator_.transfer_term(b);
    Term int_itp;
    Result r = interpolator_->get_interpolant(
        interpolator_->make_term(And, int_b, int_tr), int_f, int_itp);
    if (!r.is_unsat()) {
      throw PonoException(
          "DAR: pairwise strengthening failed, "
          "expect UNSAT interpolation query");
    }
    Term itp = to_solver_.transfer_term(int_itp);
    logger.log(3,
               "DAR: strengthening backward reachability sequence at position "
               "{} with {}",
               i + 1,
               itp);
    backward_seq_.at(i + 1) =
        solver_->make_term(And, backward_seq_.at(i + 1), unroller_.untime(itp));
  }
}

bool DualApproxReach::check_entail(const Term & p, const Term & q)
{
  solver_->reset_assertions();
  solver_->assert_formula(
      solver_->make_term(And, p, solver_->make_term(Not, q)));
  Result r = solver_->check_sat();
  assert(r.is_unsat() || r.is_sat());
  return r.is_unsat();
}

bool DualApproxReach::check_fixed_point(const TermVec & reach_seq)
{
  assert(reach_seq.size() > 1);
  Term acc_img = reach_seq.at(0);
  for (size_t i = 1; i < reach_seq.size(); ++i) {
    if (check_entail(reach_seq.at(i), acc_img)) {
      invar_ = acc_img;
      return true;
    }
    acc_img = solver_->make_term(Or, acc_img, reach_seq.at(i));
  }
  return false;
}

}  // namespace pono
