/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann
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

#include "bmc_simplepath.h"

#include "utils/logger.h"

using namespace smt;

namespace pono {

BmcSimplePath::BmcSimplePath(const Property & p, SolverEnum se) : super(p, se)
{
}

BmcSimplePath::BmcSimplePath(const Property & p, SmtSolver solver)
    : super(p, solver)
{
}

BmcSimplePath::BmcSimplePath(const PonoOptions & opt,
                             const Property & p,
                             SolverEnum se)
  : super(opt, p, se)
{
}

BmcSimplePath::BmcSimplePath(const PonoOptions & opt,
                             const Property & p,
                             smt::SmtSolver solver)
    : super(opt, p, solver)
{
}

BmcSimplePath::~BmcSimplePath() {}

ProverResult BmcSimplePath::check_until(int k)
{
  for (int i = 0; i <= k; ++i) {
    logger.log(1, "Checking Bmc at bound: {}", i);
    if (!base_step(i)) {
      return ProverResult::FALSE;
    }
    logger.log(1, "Checking simple path at bound: {}", i);
    if (cover_step(i)) {
      return ProverResult::TRUE;
    }
  }
  return ProverResult::UNKNOWN;
}

bool BmcSimplePath::cover_step(int i)
{
  if (i <= reached_k_) {
    return false;
  }

  solver_->push();
  solver_->assert_formula(init0_);
  Term not_init = solver_->make_term(PrimOp::Not, ts_.init());
  for (int j = 1; j <= i; ++j) {
    solver_->assert_formula(unroller_.at_time(not_init, j));
  }
  if (ts_.statevars().size() && check_simple_path_lazy(i)) {
    return true;
  }
  solver_->pop();

  ++reached_k_;

  return false;
}

}  // namespace pono
