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

bool DualApproxReach::step(int i) { throw PonoException("NYI"); }

bool DualApproxReach::step_0() { throw PonoException("NYI"); }

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
