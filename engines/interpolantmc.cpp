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

#include "smt-switch/exceptions.h"

#include "available_solvers.h"
#include "interpolantmc.h"

#include "utils/logger.h"

using namespace smt;

namespace pono {

InterpolantMC::InterpolantMC(Property & p, SolverEnum se)
    : super(p, se),
      interpolator_(create_interpolating_solver(se)),
      to_interpolator_(interpolator_),
      to_solver_(solver_)
{
  initialize();
}

InterpolantMC::InterpolantMC(Property & p,
                             const SmtSolver & slv,
                             const SmtSolver & itp)
    : super(p, slv),
      interpolator_(itp),
      to_interpolator_(interpolator_),
      to_solver_(solver_)
{
  initialize();
}

InterpolantMC::InterpolantMC(const PonoOptions & opt,
                             Property & p,
                             SolverEnum se)
  : super(opt, p, se),
    interpolator_(create_interpolating_solver(se)),
    to_interpolator_(interpolator_),
    to_solver_(solver_)
{
  initialize();
}

InterpolantMC::InterpolantMC(const PonoOptions & opt,
                             Property & p,
                             const SmtSolver & slv,
                             const SmtSolver & itp)
    : super(opt, p, slv),
      interpolator_(itp),
      to_interpolator_(itp),
      to_solver_(slv)
{
  initialize();
}

InterpolantMC::~InterpolantMC() {}

void InterpolantMC::initialize()
{
  super::initialize();

  reset_assertions(interpolator_);

  // symbols are already created in solver
  // need to add symbols at time 1 to cache
  // (only time 1 because Craig Interpolant has to share symbols between A and
  // B)
  UnorderedTermMap & cache = to_solver_.get_cache();
  Term tmp1;
  for (auto s : ts_.statevars()) {
    tmp1 = unroller_.at_time(s, 1);
    cache[to_interpolator_.transfer_term(tmp1)] = tmp1;
  }
  for (auto i : ts_.inputvars()) {
    tmp1 = unroller_.at_time(i, 1);
    cache[to_interpolator_.transfer_term(tmp1)] = tmp1;
  }

  concrete_cex_ = false;
  init0_ = unroller_.at_time(ts_.init(), 0);
  transA_ = unroller_.at_time(ts_.trans(), 0);
  transB_ = solver_->make_term(true);
}

ProverResult InterpolantMC::check_until(int k)
{
  try {
    for (int i = 0; i <= k; ++i) {
      if (step(i)) {
        return ProverResult::TRUE;
      } else if (concrete_cex_) {
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
  Term int_bad = to_interpolator_.transfer_term(bad_i);
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
        interpolator_->make_term(And, int_transB, int_bad),
        int_Ri);

    got_interpolant = r.is_unsat();

    if (got_interpolant) {
      Ri = to_solver_.transfer_term(int_Ri);
      // map Ri to time 0
      Ri = unroller_.at_time(unroller_.untime(Ri), 0);

      if (check_entail(Ri, R)) {
        // check if the over-approximation has reached a fix-point
        logger.log(1, "Found a proof at bound: {}", i);
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
      reset_assertions(solver_);

      Term solver_trans = solver_->make_term(And, transA_, transB_);
      solver_->assert_formula(solver_->make_term(
          And, init0_, solver_->make_term(And, solver_trans, bad_i)));

      Result r = solver_->check_sat();
      if (!r.is_sat()) {
        throw PonoException("Internal error: Expecting satisfiable result");
      }
      return false;
    }
    else if (r.is_unknown())
    {
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
  reset_assertions(solver_);
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

void InterpolantMC::reset_assertions(SmtSolver & s)
{
  // reset assertions is not supported by all solvers
  // but MathSAT is the only supported solver that can do interpolation
  // so this should be safe
  try {
    s->reset_assertions();
  }
  catch (NotImplementedException & e) {
    throw PonoException("Got unexpected solver in InterpolantMC.");
  }
}

bool InterpolantMC::check_entail(const Term & p, const Term & q)
{
  reset_assertions(solver_);
  solver_->assert_formula(
      solver_->make_term(And, p, solver_->make_term(Not, q)));
  Result r = solver_->check_sat();
  assert(r.is_unsat() || r.is_sat());
  return r.is_unsat();
}

}  // namespace pono
