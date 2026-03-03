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
 ** \brief Implementation of interpolation-sequence based model checking.
 **
 ** This unbounded model checking algorithm was introduced by Yakir Vizel and
 ** Orna Grumberg in FMCAD 2009 (https://doi.org/10.1109/FMCAD.2009.5351148).
 **
 **/

#include "engines/interp_seq_mc.h"

#include "smt-switch/exceptions.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/timestamp.h"

using namespace smt;

namespace pono {

InterpSeqMC::InterpSeqMC(const SafetyProperty & p,
                         const TransitionSystem & ts,
                         const SmtSolver & solver,
                         PonoOptions opt,
                         Engine engine)
    : super(p, ts, solver, opt, engine),
      interpolator_(
          create_interpolating_solver_for(options_.smt_interpolator_,
                                          engine_,
                                          options_.printing_smt_interpolator_,
                                          options_.smt_interpolator_opts_)),
      to_interpolator_(interpolator_),
      to_solver_(solver_)
{
}

void InterpSeqMC::initialize()
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
  reach_seq_.clear();
  trans_seq_.clear();
  int_trans_seq_.clear();
  reach_seq_.push_back(ts_.init());
}

void InterpSeqMC::reset_env()
{
  // (Soft-)Reset the solver here
  // (assuming no assertions were push to context level 0).
  // This step is needed if the engine is used in a CAGAR loop
  // as there might be assertions from previous iterations.
  assert(!concrete_cex_ || solver_->get_context_level() == 1);
  solver_->pop(solver_->get_context_level());
  // Reinitialize the prover
  initialized_ = false;
  InterpSeqMC::initialize();
}

ProverResult InterpSeqMC::check_until(int k)
{
  initialize();

  try {
    for (int i = 0; i <= k; ++i) {
      const bool step_result = step(i);
#ifndef NDEBUG
      const std::size_t log_level = (step_result || concrete_cex_) ? 2 : 4;
      logger.log(log_level,
                 "Interpolation stats: {} calls took {:.3f} s",
                 total_interp_call_count_,
                 total_interp_call_time_);
#endif
      if (step_result) {
        return ProverResult::TRUE;
      } else if (concrete_cex_) {
        compute_witness();
        return ProverResult::FALSE;
      }
    }
  }
  catch (InternalSolverException & e) {
    logger.log(1, "ISMC failed due to: {}", e.what());
  }
  return ProverResult::UNKNOWN;
}

bool InterpSeqMC::step(int i)
{
  if (i <= reached_k_) {
    return false;
  }

  logger.log(1, "Running ISMC at bound: {}", i);

  if (i == 0) {
    // Can't get an interpolant at bound 0
    // only checking for trivial bug
    return step_0();
  }

  // construct the transition formula at current time step
  update_term_map(i);
  const smt::Term trans_i = unroller_.at_time(ts_.trans(), i - 1);
  trans_seq_.push_back(trans_i);
  if (i == 1) {
    // for convenience, we conjoin TR(0, 1) with Init(0)
    assert(trans_seq_.size() == 1);
    trans_seq_.at(0) =
        solver_->make_term(And, unroller_.at_time(ts_.init(), 0), trans_i);
  }
  int_trans_seq_.push_back(to_interpolator_.transfer_term(trans_seq_.back()));

  Term bad_i = unroller_.at_time(bad_, i);
  Term int_bad_i = to_interpolator_.transfer_term(bad_i);

  // temporarily push `bad` to `trans_seq` to avoid copying the whole vector
  int_trans_seq_.push_back(int_bad_i);

#ifndef NDEBUG
  const std::clock_t start_t = std::clock();
#endif
  TermVec int_itp_seq;
  Result r =
      interpolator_->get_sequence_interpolants(int_trans_seq_, int_itp_seq);

#ifndef NDEBUG
  log_interp_time(start_t, total_interp_call_count_, total_interp_call_time_);
  if (r.is_unsat()) {
    check_itp_sequence(int_trans_seq_, int_itp_seq);
  }
#endif

  // pop `bad` out from `trans_seq`
  int_trans_seq_.pop_back();

  if (r.is_unsat()) {
    // update reachability sequence with interpolants
    reach_seq_.push_back(solver_->make_term(true));
    assert(reach_seq_.size() == int_itp_seq.size() + 1);
    for (size_t j = 0; j < int_itp_seq.size(); ++j) {
      Term itp = unroller_.untime(to_solver_.transfer_term(int_itp_seq.at(j)));
      reach_seq_.at(j + 1) = solver_->make_term(And, reach_seq_.at(j + 1), itp);
    }
    if (check_fixed_point()) {
      logger.log(1, "Found a proof at bound: {}", i);
      return true;
    }
  } else {  // handle sat and unknown
    // check the BMC query using the solver with model generation
    // note that we also perform the check even when interpolation fails
    // (i.e., r.is_unknown()), because iterative construction of interpolation
    // sequence may fail if the given formula is satisfiable
    assert(solver_->get_context_level() == 0);
    solver_->push();
    Term trans_until_i = (trans_seq_.size() == 1)
                             ? trans_seq_.at(0)
                             : solver_->make_term(And, trans_seq_);
    solver_->assert_formula(solver_->make_term(And, trans_until_i, bad_i));
    Result bmc_res = solver_->check_sat();
    if (bmc_res.is_sat()) {
      // found a concrete counter example
      // replay it in the solver with model generation
      concrete_cex_ = true;
    } else if (r.is_sat() && !bmc_res.is_sat()) {
      throw PonoException("Internal error: Expecting satisfiable result");
    } else {
      throw PonoException("Interpolation failed due to: "
                          + r.get_explanation());
    }
    return false;
  }

  ++reached_k_;
  return false;
}

bool InterpSeqMC::step_0()
{
  assert(solver_->get_context_level() == 0);
  // push the unrolled formulas here
  // as compute_witness() rely on timed variables
  solver_->push();
  solver_->assert_formula(unroller_.at_time(reach_seq_.at(0), 0));
  solver_->assert_formula(unroller_.at_time(bad_, 0));

  Result r = solver_->check_sat();
  if (r.is_unsat()) {
    reached_k_ = 0;
    solver_->pop();
  } else if (r.is_sat()) {
    // do not pop here to keep the solver state
    // for later witness extraction (`compute_witness()`)
    concrete_cex_ = true;
  } else {
    throw PonoException("ISMC: step_0 failed, unexpected result: "
                        + r.to_string());
  }
  return false;
}

void InterpSeqMC::update_term_map(size_t i)
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

bool InterpSeqMC::check_fixed_point()
{
  assert(reach_seq_.size() > 1);
  assert(solver_->get_context_level() == 0);
  // initialize solver stack and reached set
  Term acc_img = reach_seq_.at(0);
  solver_->push();
  solver_->assert_formula(solver_->make_term(Not, acc_img));
  for (size_t i = 1; i < reach_seq_.size(); ++i) {
    solver_->push();
    solver_->assert_formula(reach_seq_.at(i));
    Result r = solver_->check_sat();
    assert(r.is_unsat() || r.is_sat());
    solver_->pop();
    if (r.is_unsat()) {
      // reached a fixed point
      solver_->pop();  // pop all assertions
      invar_ = acc_img;
      return true;
    }
    // extend the accumulated reached set
    solver_->assert_formula(solver_->make_term(Not, reach_seq_.at(i)));
    acc_img = solver_->make_term(Or, acc_img, reach_seq_.at(i));
  }
  solver_->pop();  // pop all assertions
  return false;
}

// check if the interpolation sequence is valid
// i.e., `A => Itp@i` holds and `Itp@i & B` is UNSAT for each i
void InterpSeqMC::check_itp_sequence(const TermVec & int_formulas,
                                     const TermVec & int_itp_seq)
{
  assert(int_formulas.size() == int_itp_seq.size() + 1);
  assert(solver_->get_context_level() == 0);
  for (size_t i = 0; i < int_itp_seq.size(); ++i) {
    TermVec int_a_vec(int_formulas.begin(), int_formulas.begin() + i + 1);
    TermVec int_b_vec(int_formulas.begin() + i + 1, int_formulas.end());
    Term int_a = (int_a_vec.size() == 1)
                     ? int_a_vec.at(0)
                     : interpolator_->make_term(And, int_a_vec);
    Term int_b = (int_b_vec.size() == 1)
                     ? int_b_vec.at(0)
                     : interpolator_->make_term(And, int_b_vec);
    Term a = to_solver_.transfer_term(int_a);
    Term b = to_solver_.transfer_term(int_b);
    Term itp = to_solver_.transfer_term(int_itp_seq.at(i));
    solver_->push();
    solver_->assert_formula(a);
    solver_->assert_formula(solver_->make_term(Not, itp));
    if (!solver_->check_sat().is_unsat()) {
      throw PonoException("Invalid interpolation sequence: A =\\> Itp@"
                          + std::to_string(i));
    }
    solver_->pop();
    solver_->push();
    solver_->assert_formula(solver_->make_term(And, itp, b));
    if (!solver_->check_sat().is_unsat()) {
      throw PonoException("Invalid interpolation sequence: Itp@"
                          + std::to_string(i) + " & B is not UNSAT");
    }
    solver_->pop();
  }
}

}  // namespace pono
