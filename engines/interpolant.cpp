#include "smt-switch/exceptions.h"

#include "interpolant.h"
#include "utils/logger.h"

using namespace smt;

namespace cosa {

InterpolantMC::InterpolantMC(const Property & p,
                             SmtSolver & slv,
                             SmtSolver & itp)
    : super(p, slv), interpolator_(itp), to_interpolator_(itp), to_solver_(slv)
{
  // symbols are already created in solver
  // need to add symbols at time 1 to cache
  // (only time 1 because Craig Interpolant has to share symbols between A and
  // B)
  UnorderedTermMap & cache = to_solver_.get_cache();
  Term tmp1;
  for (auto s : ts_.states())
  {
    unroller_.at_time(s, 0);
    tmp1 = unroller_.at_time(s, 1);
    cache[to_interpolator_.transfer_term(tmp1)] = tmp1;
  }
  for (auto i : ts_.inputs())
  {
    unroller_.at_time(i, 0);
    tmp1 = unroller_.at_time(i, 1);
    cache[to_interpolator_.transfer_term(tmp1)] = tmp1;
  }
  initialize();
}

InterpolantMC::~InterpolantMC() {}

void InterpolantMC::initialize()
{
  super::initialize();

  concrete_cex_ = false;

  // reset assertions is not supported by all solvers
  // but MathSAT is the only supported solver that can do interpolation
  // so this should be safe
  try
  {
    interpolator_->reset_assertions();
  }
  catch (NotImplementedException & e)
  {
    throw CosaException("Got unexpected solver in InterpolantMC.");
  }

  // populate map from time 1 to time 0
  for (auto s : ts_.states())
  {
    Term s0 = unroller_.at_time(s, 0);
  }

  for (auto i : ts_.inputs())
  {
    Term i0 = unroller_.at_time(i, 0);
  }

  init0_ = unroller_.at_time(ts_.init(), 0);
  R_ = init0_;
  Ri_ = solver_->make_term(true);
  transA_ = unroller_.at_time(ts_.trans(), 0);
  transB_ = solver_->make_term(true);
}

ProverResult InterpolantMC::check_until(int k)
{
  try
  {
    for (int i = 0; i <= k; ++i)
    {
      if (step(i))
      {
        return ProverResult::TRUE;
      }
      else if (concrete_cex_)
      {
        return ProverResult::FALSE;
      }
    }
  }
  catch (InternalSolverException & e)
  {
    logger.log(1, "Failed when computing interpolant.");
  }
  return ProverResult::UNKNOWN;
}

bool InterpolantMC::step(int i)
{
  if (i <= reached_k_)
  {
    return false;
  }

  logger.log(1, "Checking interpolation at bound: {}", i);
  Term bad_i = unroller_.at_time(bad_, i);

  // These bools are always opposite except at i = 0
  bool got_interpolant = true;
  bool is_sat = false;

  R_ = init0_;

  while (got_interpolant)
  {
    if (i > 0)
    {
      Term int_R = to_interpolator_.transfer_term(R_);
      Term int_transA = to_interpolator_.transfer_term(transA_);
      Term int_transB = to_interpolator_.transfer_term(transB_);
      Term int_bad = to_interpolator_.transfer_term(bad_i);
      Term int_Ri;
      got_interpolant = interpolator_->get_interpolant(
          interpolator_->make_term(And, int_R, int_transA),
          interpolator_->make_term(And, int_transB, int_bad),
          int_Ri);

      if (got_interpolant)
      {
        Ri_ = to_solver_.transfer_term(int_Ri);
      }

      is_sat = !got_interpolant;
    }
    else
    {
      // Can't get an interpolant at bound 0
      // only checking for trivial bug
      solver_->reset_assertions();

      solver_->push();
      solver_->assert_formula(R_);
      solver_->assert_formula(bad_i);
      Result r = solver_->check_sat();
      solver_->pop();

      got_interpolant = false;
      is_sat = r.is_sat();
    }

    if (is_sat && (R_ == init0_))
    {
      // found a concrete counter example
      // replay it in the solver with model generation
      concrete_cex_ = true;
      solver_->reset_assertions();
      Term solver_trans = solver_->make_term(And, transA_, transB_);
      solver_->assert_formula(solver_->make_term(
          And, init0_, solver_->make_term(And, solver_trans, bad_i)));
      Result r = solver_->check_sat();
      if (!r.is_sat())
      {
        throw CosaException("Internal error: Expecting satisfiable result");
      }
      ++reached_k_;
      return false;
    }
    else if (got_interpolant)
    {
      // map Ri to time 0
      Ri_ = unroller_.at_time(unroller_.untime(Ri_), 0);

      if (check_overapprox())
      {
        logger.log(1, "Found a proof at bound: {}", i);
        return true;
      }
      else
      {
        logger.log(1, "Extending initial states.");
        R_ = solver_->make_term(Or, R_, Ri_);
      }
    }
  }

  // Note: important that it's for i > 0
  // transB can't have any symbols from time 0 in it
  if (i > 0)
  {
    // extend the unrolling
    transB_ =
        solver_->make_term(And, transB_, unroller_.at_time(ts_.trans(), i));
  }

  ++reached_k_;
  return false;
}

bool InterpolantMC::check_overapprox()
{
  solver_->reset_assertions();
  Term Rp = R_;
  Term Rpi = Ri_;
  solver_->assert_formula(
      solver_->make_term(And, Rpi, solver_->make_term(Not, Rp)));
  Result r = solver_->check_sat();
  if (r.is_unsat())
  {
    return true;
  }
  else
  {
    return false;
  }
}

}  // namespace cosa
