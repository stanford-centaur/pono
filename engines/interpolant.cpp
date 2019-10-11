#include "smt-switch/exceptions.h"

#include "interpolant.h"
#include "utils/logger.h"

using namespace smt;

namespace cosa {

InterpolantMC::InterpolantMC(const Property & p,
                             smt::SmtSolver & itp,
                             smt::SmtSolver & slv)
    : ts_(p.transition_system()),
      property_(p),
      interpolator_(itp),
      solver_(slv),
      unroller_(ts_, interpolator_)
{
  initialize();
}

InterpolantMC::~InterpolantMC() {}

void InterpolantMC::initialize()
{
  reached_k_ = -1;
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
    map_1_to_0[unroller_.at_time(s, 1)] = s0;
  }

  for (auto i : ts_.inputs())
  {
    Term i0 = unroller_.at_time(i, 0);
    map_1_to_0[unroller_.at_time(i, 1)] = i0;
  }

  bad_ = interpolator_->make_term(Not, property_.prop());
  init0_ = unroller_.at_time(ts_.init(), 0);
  R_ = init0_;
  Ri_ = interpolator_->make_value(true);
  transA_ = unroller_.at_time(ts_.trans(), 0);
  transB_ = interpolator_->make_value(true);
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

ProverResult InterpolantMC::prove()
{
  try
  {
    for (int i = 0;; ++i)
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

bool InterpolantMC::witness(std::vector<UnorderedTermMap> & out)
{
  // TODO: make sure the solver state is SAT

  for (int i = 0; i <= reached_k_; ++i)
  {
    out.push_back(UnorderedTermMap());
    UnorderedTermMap & map = out.back();

    for (auto v : ts_.states())
    {
      Term vi = solver_->transfer_term(unroller_.at_time(v, i));
      Term r = solver_->get_value(vi);
      map[v] = r;
    }

    for (auto v : ts_.inputs())
    {
      Term vi = solver_->transfer_term(unroller_.at_time(v, i));
      Term r = solver_->get_value(vi);
      map[v] = r;
    }
  }

  return true;
}

bool InterpolantMC::step(int i)
{
  if (i <= reached_k_)
  {
    return false;
  }

  logger.log(1, "Checking interpolation at bound: {}", i);

  Term bad_i = unroller_.at_time(bad_, i);
  bool got_interpolant = true;

  R_ = init0_;

  while (got_interpolant)
  {
    if (i > 0)
    {
      got_interpolant = interpolator_->get_interpolant(
          interpolator_->make_term(And, R_, transA_),
          interpolator_->make_term(And, transB_, bad_i),
          Ri_);
    }
    else
    {
      got_interpolant = interpolator_->get_interpolant(R_, bad_i, Ri_);
    }

    if (!got_interpolant && (R_ == init0_))
    {
      // found a concrete counter example
      // replay it in the solver with model generation
      concrete_cex_ = true;
      Term solver_init = solver_->transfer_term(init0_);
      Term solver_trans = solver_->make_term(And,
                                             solver_->transfer_term(transA_),
                                             solver_->transfer_term(transB_));
      Term solver_bad = solver_->transfer_term(bad_i);
      solver_->assert_formula(solver_->make_term(
          And, solver_init, solver_->make_term(And, solver_trans, solver_bad)));
      Result r = solver_->check_sat();
      if (!r.is_sat())
      {
        throw CosaException("Internal error: Expecting satisfiable result");
      }
      return false;
    }
    else if (got_interpolant)
    {
      // map Ri to time 0
      Ri_ = interpolator_->substitute(Ri_, map_1_to_0);

      if (check_overapprox())
      {
        logger.log(1, "Found a proof at bound: {}", i);
        return true;
      }
      else
      {
        logger.log(1, "Extending initial states.");
        R_ = interpolator_->make_term(Or, R_, Ri_);
      }
    }
  }

  // Note: important that it's for i > 0
  // transB can't have any symbols from time 0 in it
  if (i > 0)
  {
    // extend the unrolling
    transB_ = interpolator_->make_term(
        And, transB_, unroller_.at_time(ts_.trans(), i));
  }

  return false;
}

bool InterpolantMC::check_overapprox()
{
  solver_->reset_assertions();
  Term Rp = solver_->transfer_term(R_);
  Term Rpi = solver_->transfer_term(Ri_);
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
