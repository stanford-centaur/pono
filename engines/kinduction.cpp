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

KInduction::KInduction(const Property & p, const TransitionSystem & ts,
                       const SmtSolver & solver,
                       PonoOptions opt)
  : super(p, ts, solver, opt)
{
  engine_ = Engine::KIND;
}

KInduction::~KInduction() {}

void KInduction::initialize()
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
  init0_ = unroller_.at_time(ts_.init(), 0);
  false_ = solver_->make_term(false);
  simple_path_ = solver_->make_term(true);
}

ProverResult KInduction::check_until(int k)
{
  initialize();
  assert(reached_k_ == -1);

  Result res;
  for (int i = reached_k_ + 1; i <= k; ++i) {

    solver_->push();

    // inductive case check
    if (!options_.kind_no_simple_path_check_)
      solver_->assert_formula(simple_path_);
    solver_->assert_formula(unroller_.at_time(bad_, i));
    logger.log(1, "Checking k-induction inductive step at bound: {}", i);
    // solver call inside 'check_simple_path_lazy/eager'
    if (!options_.kind_eager_simple_path_check_) {
      if (ts_.statevars().size() && check_simple_path_lazy(i)) {
	return ProverResult::TRUE;
      }
    } else {
      if (ts_.statevars().size() && check_simple_path_eager(i)) {
	return ProverResult::TRUE;
      }
    }

//DELETE THIS
//    res = solver_->check_sat();
//    if (res.is_unsat()) {
//      return ProverResult::TRUE;
//    }

    // base case check
    solver_->assert_formula(init0_);
    logger.log(1, "Checking k-induction base case at bound: {}", i);
    res = solver_->check_sat();
    if (res.is_sat()) {
      compute_witness();
      return ProverResult::FALSE;
    }

    solver_->pop();

    // add transition and negated bad state property
    solver_->assert_formula(unroller_.at_time(ts_.trans(), i));
    solver_->assert_formula(unroller_.at_time(solver_->make_term(Not, bad_), i));

    reached_k_++;
    
#if 0    
    ///////////////////////////////////
    logger.log(1, "Checking k-induction base case at bound: {}", i);
    if (!base_step(i)) {
      compute_witness();
      return ProverResult::FALSE;
    }
    logger.log(1, "Checking k-induction inductive step at bound: {}", i);
    if (inductive_step(i)) {
      return ProverResult::TRUE;
    }
    /////////////////////////////////// 
 #endif
    
  } //end: for all bounds
  
  return ProverResult::UNKNOWN;
}

bool KInduction::base_step(int i)
{
  //TODO function deprecated
  abort();

  if (i <= reached_k_) {
    return true;
  }

  solver_->push();
  solver_->assert_formula(init0_);
  solver_->assert_formula(unroller_.at_time(bad_, i));
  Result r = solver_->check_sat();
  if (r.is_sat()) {
    return false;
  }
  solver_->pop();

  const Term &prop = solver_->make_term(Not, bad_);
  solver_->assert_formula(unroller_.at_time(ts_.trans(), i));
  solver_->assert_formula(unroller_.at_time(prop, i));

  return true;
}

bool KInduction::inductive_step(int i)
{
  //TODO function deprecated
  abort();

  if (i <= reached_k_) {
    return false;
  }

  solver_->push();
  if (!options_.kind_no_simple_path_check_)
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
  assert(!options_.kind_no_simple_path_check_);
  assert(ts_.statevars().size());

  Term disj = false_;
  for (const auto &v : ts_.statevars()) {
    Term vi = unroller_.at_time(v, i);
    Term vj = unroller_.at_time(v, j);
    Term eq = solver_->make_term(PrimOp::Equal, vi, vj);
    Term neq = solver_->make_term(PrimOp::Not, eq);
    disj = solver_->make_term(PrimOp::Or, disj, neq);
  }

  return disj;
}

bool KInduction::check_simple_path_eager(int i)
{
  assert(options_.kind_eager_simple_path_check_);
  logger.log(2, "  Eagerly checking k-induction simple path at bound: {}", i);

  const bool no_simp_path_check = options_.kind_no_simple_path_check_;

  //Here we assume that all pairs for values smaller than 'i' have been added
  //If simple path checking is disabled then we still need the final
  //solver call below for inductive case check
  for (int j = 0; (!no_simp_path_check && j < i); j++) {
    Term constraint = simple_path_constraint(j, i);
    logger.log(3, "   Adding simple path clause for pair 'j,i' = {},{}", j,i);
    simple_path_ =
      solver_->make_term(PrimOp::And, simple_path_, constraint);
    solver_->assert_formula(constraint);
  }

  logger.log(2, "    Calling solver for simple path check");
  Result r = solver_->check_sat();
  if (r.is_unsat()) {
    logger.log(2, "      Simple path check UNSAT");
    return true;
  }

  return false;
}

bool KInduction::check_simple_path_lazy(int i)
{
  logger.log(2, "  Lazily checking k-induction simple path at bound: {}", i);
  bool added_to_simple_path = false;

  do {
    logger.log(2, "    Calling solver for simple path check");
    Result r = solver_->check_sat();
    if (r.is_unsat()) {
      logger.log(2, "      Simple path check UNSAT");
      return true;
    }

    // We need the above check_sat call if simple path checks are
    // disabled because it is part of the inductive step check. Hence
    // we return immediately here and skip simple path checking.
    if (options_.kind_no_simple_path_check_)
      return false;

    added_to_simple_path = false;

    for (int j = 0; j < i && !added_to_simple_path; ++j) {
      for (int l = j + 1; l <= i; ++l) {
        Term constraint = simple_path_constraint(j, l);
	logger.log(3, "    Checking constraint for pair j,l = {} , {}", j,l);
        if (solver_->get_value(constraint) == false_) {
	  logger.log(3, "      Adding constraint for pair j,l = {} , {}", j,l);
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
