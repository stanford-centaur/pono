/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#include "engines/interpolantmc.h"

#include "smt-switch/exceptions.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"

using namespace smt;

namespace pono {

InterpolantMC::InterpolantMC(const Property & p,
                             const TransitionSystem & ts,
                             const SmtSolver & slv,
                             PonoOptions opt)
    : super(p, ts, slv, opt),
      interpolator_(create_interpolating_solver_for(opt.smt_interpolator_,
                                                    Engine::INTERP)),
      to_interpolator_(interpolator_),
      to_solver_(solver_)
{
  engine_ = Engine::INTERP;
}

InterpolantMC::~InterpolantMC() {}

void InterpolantMC::initialize()
{
  if (initialized_) {
    return;
  }

  super::initialize();

  interpolator_->reset_assertions();

  // symbols are already created in solver
  // need to add symbols at time 1 to cache
  // (only at time 1 because Craig Interpolation has to share symbols
  // between A and B)
  UnorderedTermMap & cache = to_solver_.get_cache();
  Term tmp1;
  for (const auto & s : ts_.statevars()) {
    tmp1 = unroller_.at_time(s, 1);
    cache[to_interpolator_.transfer_term(tmp1)] = tmp1;
  }
  for (const auto & i : ts_.inputvars()) {
    tmp1 = unroller_.at_time(i, 1);
    cache[to_interpolator_.transfer_term(tmp1)] = tmp1;
  }

  // need to copy over UF as well
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
  init0_ = unroller_.at_time(ts_.init(), 0);
  transA_ = unroller_.at_time(ts_.trans(), 0);
  transB_ = solver_->make_term(true);
  bad_disjuncts_ = solver_->make_term(false);
}

ProverResult InterpolantMC::check_until(int k)
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
    logger.log(1, "Failed when computing interpolant.");
  }
  return ProverResult::UNKNOWN;
}

bool InterpolantMC::step(int i)
{
  if (i <= reached_k_) {
    return false;
  }

  logger.log(1, "Checking interpolation at bound: {}", i);

  if (i == 0) {
    // Can't get an interpolant at bound 0
    // only checking for trivial bug
    return step_0();
  }

  Term bad_i = unroller_.at_time(bad_, i);
  bad_disjuncts_ = solver_->make_term(Or, bad_disjuncts_, bad_i);
  Term int_bad_disjuncts = to_interpolator_.transfer_term(bad_disjuncts_);
  Term int_transA = to_interpolator_.transfer_term(transA_);
  Term int_transB = to_interpolator_.transfer_term(transB_);
  Term R = init0_;
  Term Ri;
  bool got_interpolant = true;

  while (got_interpolant) {
    Term int_R = to_interpolator_.transfer_term(R);
    Term int_Ri;
    Result r = interpolator_->get_interpolant(
        interpolator_->make_term(And, int_R, int_transA),
        interpolator_->make_term(And, int_transB, int_bad_disjuncts),
        int_Ri);

    got_interpolant = r.is_unsat();

    if (got_interpolant) {
      Ri = to_solver_.transfer_term(int_Ri);
      // map Ri to time 0
      Ri = unroller_.at_time(unroller_.untime(Ri), 0);

      if (check_entail(Ri, R)) {
        // check if the over-approximation has reached a fix-point
        logger.log(1, "Found a proof at bound: {}", i);
        invar_ = unroller_.untime(R);
        return true;
      } else {
        logger.log(1, "Extending initial states.");
        logger.log(3, "Using interpolant: {}", Ri);
        R = solver_->make_term(Or, R, Ri);
      }
    } else if (R == init0_) {
      // found a concrete counter example
      // replay it in the solver with model generation
      concrete_cex_ = true;
      solver_->reset_assertions();

      Term solver_trans = solver_->make_term(And, transA_, transB_);
      solver_->assert_formula(solver_->make_term(
          And, init0_, solver_->make_term(And, solver_trans, bad_i)));

      Result r = solver_->check_sat();
      if (!r.is_sat()) {
        throw PonoException("Internal error: Expecting satisfiable result");
      }
      return false;
    } else if (r.is_unknown()) {
      // TODO: figure out if makes sense to increase bound and try again
      throw PonoException("Interpolant generation failed.");
    }
  }

  // Note: important that it's for i > 0
  // transB can't have any symbols from time 0 in it
  assert(i > 0);
  // extend the unrolling
  transB_ = solver_->make_term(And, transB_, unroller_.at_time(ts_.trans(), i));
  ++reached_k_;

  return false;
}

bool InterpolantMC::step_0()
{
  solver_->reset_assertions();
  solver_->assert_formula(init0_);
  solver_->assert_formula(unroller_.at_time(bad_, 0));

  Result r = solver_->check_sat();
  if (r.is_unsat()) {
    ++reached_k_;
  } else {
    concrete_cex_ = true;
  }
  return false;
}

bool InterpolantMC::check_entail(const Term & p, const Term & q)
{
  solver_->reset_assertions();
  solver_->assert_formula(
      solver_->make_term(And, p, solver_->make_term(Not, q)));
  Result r = solver_->check_sat();
  assert(r.is_unsat() || r.is_sat());
  return r.is_unsat();
}

}  // namespace pono
