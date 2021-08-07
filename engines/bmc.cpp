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

#include "bmc.h"
#include "utils/logger.h"

using namespace smt;

namespace pono {

Bmc::Bmc(const Property & p, const TransitionSystem & ts,
         const SmtSolver & solver, PonoOptions opt)
  : super(p, ts, solver, opt)
{
  engine_ = Engine::BMC;
}

Bmc::~Bmc() {}

void Bmc::initialize()
{
  if (initialized_) {
    return;
  }

  super::initialize();

  // NOTE: There's an implicit assumption that this solver is only used for
  // model checking once Otherwise there could be conflicting assertions to
  // the solver or it could just be polluted with redundant assertions in the
  // future we can use solver_->reset_assertions(), but it is not currently
  // supported in boolector
  solver_->assert_formula(unroller_.at_time(ts_.init(), 0));
}

ProverResult Bmc::check_until(int k)
{
  initialize();

  const int step_bound = 1;
  const int start_bound = 0;

  // reached_k == -1 initially
  
  //  if (start_bound > 0)
  //  reached_k_ = start_bound - 1;
  
//  for (int j = 1; j < start_bound; ++j)
  // {
  //   std::cout << "DEBUG (check-until) adding trans for j-1 == " << j - 1 << std::endl;
  //   solver_->assert_formula(unroller_.at_time(ts_.trans(), j - 1));
//  }
  
  for (int i = start_bound; i <= k; i+=step_bound /*i = i == 0 ? 1 : i << 1*/) {
    if (!step(i)) {
      compute_witness();
      return ProverResult::FALSE;
    }
  }
  return ProverResult::UNKNOWN;
}

bool Bmc::step(int i)
{
  if (i <= reached_k_) {
    return true;
  }

  bool res = true;
  if (i > 0) {
//
    std::cout << "DEBUG reached k " << reached_k_ << ", i " << i << std::endl;
    for (int j = reached_k_ == -1 ? 1 : reached_k_ + 1; j <= i; j++) {
      std::cout << "DEBUG adding trans for j-1 == " << j - 1 << std::endl;
      solver_->assert_formula(unroller_.at_time(ts_.trans(), j - 1));
    }
    //OLD
//    solver_->assert_formula(unroller_.at_time(ts_.trans(), i - 1));
  }

  solver_->push();
  logger.log(1, "Checking bmc at bound: {}", i);

  const int cex_guarantee = 1;

  Term clause;
  if (cex_guarantee) {
    // make sure we cover *all* states
  // TODO (not critical): can make 'solver_->make_term(false)' a constant in the object
    clause = solver_->make_term(false);
    for (int j = reached_k_ + 1; j <= i; j++) {
      std::cout << "DEBUG adding bad state literal for j == " << j << std::endl;
      clause = solver_->make_term(PrimOp::Or, clause, unroller_.at_time(bad_, j));
    }
  } else {
    clause = unroller_.at_time(bad_, i);
  }
  
  solver_->assert_formula(clause);
  
  //TODO: add bad state clause here
    // OLD
    // solver_->assert_formula(unroller_.at_time(bad_, i));


  Result r = solver_->check_sat();
  if (r.is_sat()) {
    res = false;
  } else {
    solver_->pop();
    reached_k_ = i;
  }

  return res;
}

}  // namespace pono
