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

#include "kinduction.h"
#include "utils/logger.h"

using namespace smt;

namespace pono {

KInduction::KInduction(const Property & p, SolverEnum se) : super(p, se)
{
  initialize();
}

KInduction::KInduction(const Property & p, const SmtSolver & solver)
    : super(p, solver)
{
  initialize();
}

KInduction::KInduction(const PonoOptions & opt,
                       const Property & p,
                       SolverEnum se)
    : super(opt, p, se)
{
  initialize();
}

KInduction::KInduction(const PonoOptions & opt,
                       const Property & p,
                       const smt::SmtSolver & solver)
    : super(opt, p, solver)
{
  initialize();
}

KInduction::~KInduction() {}

void KInduction::initialize()
{
  super::initialize();
  // NOTE: There's an implicit assumption that this solver is only used for
  // model checking once Otherwise there could be conflicting assertions to
  // the solver or it could just be polluted with redundant assertions in the
  // future we can use solver_->reset_assertions(), but it is not currently
  // supported in boolector
  init0_ = unroller_.at_time(ts_.init(), 0);
  false_ = solver_->make_term(false);
  simple_path_ = solver_->make_term(true);
}

ProverResult KInduction::check_until(int k)
{
  for (int i = 0; i <= k; ++i) {
    logger.log(1, "Checking k-induction base case at bound: {}", i);
    if (!base_step(i)) {
      return ProverResult::FALSE;
    }
    logger.log(1, "Checking k-induction inductive step at bound: {}", i);
    if (inductive_step(i)) {
      return ProverResult::TRUE;
    }
  }
  return ProverResult::UNKNOWN;
}

bool KInduction::base_step(int i)
{
  if (i <= reached_k_) {
    return true;
  }

  solver_->push();
  solver_->assert_formula(init0_);
  solver_->assert_formula(unroller_.at_time(bad_, i));
  Result r = solver_->check_sat();
  if (r.is_sat()) {
    ++reached_k_;
    return false;
  }
  solver_->pop();

  solver_->assert_formula(unroller_.at_time(ts_.trans(), i));
  solver_->assert_formula(unroller_.at_time(property_.prop(), i));

  return true;
}

bool KInduction::inductive_step(int i)
{
  if (i <= reached_k_) {
    return false;
  }

  solver_->push();
  solver_->assert_formula(simple_path_);
  solver_->assert_formula(unroller_.at_time(bad_, i + 1));

  if (ts_.statevars().size() && check_simple_path_lazy(i + 1)) {
    return true;
  }

  solver_->pop();

  ++reached_k_;

  return false;
}

Term KInduction::simple_path_constraint(int i, int j)
{
  assert(ts_.statevars().size());

  Term disj = false_;
  for (auto v : ts_.statevars()) {
    Term vi = unroller_.at_time(v, i);
    Term vj = unroller_.at_time(v, j);
    Term eq = solver_->make_term(PrimOp::Equal, vi, vj);
    Term neq = solver_->make_term(PrimOp::Not, eq);
    disj = solver_->make_term(PrimOp::Or, disj, neq);
  }

  return disj;
}

bool KInduction::check_simple_path_lazy(int i)
{
  bool added_to_simple_path = false;

  do {
    Result r = solver_->check_sat();
    if (r.is_unsat()) {
      return true;
    }

    added_to_simple_path = false;

    for (int j = 0; j < i && !added_to_simple_path; ++j) {
      for (int l = j + 1; l <= i; ++l) {
        Term constraint = simple_path_constraint(j, l);
        if (solver_->get_value(constraint) == false_) {
	  logger.log(2, "Adding Simple Path Clause");
          simple_path_ =
              solver_->make_term(PrimOp::And, simple_path_, constraint);
          solver_->assert_formula(constraint);
          added_to_simple_path = true;
          break;
        }
      }
    }

  } while (added_to_simple_path);

  return false;
}

}  // namespace pono
