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
#include "utils/term_analysis.h"

using namespace smt;

namespace pono {

InterpSeqMC::InterpSeqMC(const SafetyProperty & p,
                         const TransitionSystem & ts,
                         const SmtSolver & slv,
                         PonoOptions opt)
    : super(p, ts, slv, opt),
      interpolator_(
          create_interpolating_solver_for(opt.smt_interpolator_, Engine::ISMC)),
      to_interpolator_(interpolator_),
      to_solver_(solver_)
{
  engine_ = Engine::ISMC;
}

InterpSeqMC::~InterpSeqMC() {}

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
  reach_seq_.push_back(ts_.init());
}

ProverResult InterpSeqMC::check_until(int k)
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

  TermVec int_itp_seq;
  Result r =
      interpolator_->get_sequence_interpolants(int_trans_seq_, int_itp_seq);

#ifndef NDEBUG
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
    solver_->reset_assertions();
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
  solver_->reset_assertions();
  // push the unrolled formulas here
  // as compute_witness() rely on timed variables
  solver_->assert_formula(unroller_.at_time(reach_seq_.at(0), 0));
  solver_->assert_formula(unroller_.at_time(bad_, 0));

  Result r = solver_->check_sat();
  if (r.is_unsat()) {
    ++reached_k_;
  } else {
    concrete_cex_ = true;
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

bool InterpSeqMC::check_entail(const Term & p, const Term & q)
{
  solver_->reset_assertions();
  solver_->assert_formula(
      solver_->make_term(And, p, solver_->make_term(Not, q)));
  Result r = solver_->check_sat();
  assert(r.is_unsat() || r.is_sat());
  return r.is_unsat();
}

bool InterpSeqMC::check_fixed_point()
{
  assert(reach_seq_.size() > 1);
  Term acc_img = reach_seq_.at(0);
  for (size_t i = 1; i < reach_seq_.size(); ++i) {
    if (check_entail(reach_seq_.at(i), acc_img)) {
      invar_ = acc_img;
      return true;
    }
    acc_img = solver_->make_term(Or, acc_img, reach_seq_.at(i));
  }
  return false;
}

// check if the interpolation sequence is valid
// i.e., `A => Itp@i` holds and `Itp@i & B` is UNSAT for each i
void InterpSeqMC::check_itp_sequence(const TermVec & int_formulas,
                                     const TermVec & int_itp_seq)
{
  assert(int_formulas.size() == int_itp_seq.size() + 1);
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
    if (!check_entail(a, itp)) {
      throw PonoException("Invalid interpolation sequence: A =\\> Itp@"
                          + std::to_string(i));
    }
    solver_->reset_assertions();
    solver_->assert_formula(solver_->make_term(And, itp, b));
    if (!solver_->check_sat().is_unsat()) {
      throw PonoException("Invalid interpolation sequence: Itp@"
                          + std::to_string(i) + " & B is not UNSAT");
    }
  }
}

}  // namespace pono
