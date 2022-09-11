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
  kind_engine_name_ = "k-induction";
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
  neg_init_terms_ = solver_->make_term(true);
  // selector literal to toggle initial state predicate
  Sort boolsort = solver_->make_sort(smt::BOOL);
  sel_init_ = solver_->make_symbol("sel_init", boolsort);
  not_sel_init_ = solver_->make_term(Not, sel_init_);
  //permanently add term '(sel_init0_ OR init0_)'
  solver_->assert_formula(solver_->make_term(PrimOp::Or, sel_init_, init0_));
}

ProverResult KInduction::check_until(int k)
{
  initialize();
  assert(reached_k_ == -1);

  Result res;
  for (int i = reached_k_ + 1; i <= k; ++i) {

    logger.log(1, "");
    kind_log_msg(1, "", "current unrolling depth/bound: {}", i);

    //disable initial state predicate
    if (!sel_assumption_.empty())
      sel_assumption_.pop_back();
    assert(sel_assumption_.size() == 0);
    sel_assumption_.push_back(sel_init_);

    // simple path check
    if (!options_.kind_no_simple_path_check_) {
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
    }

    if (i >= 1 && !options_.kind_no_ind_check_init_states_ &&
       !options_.kind_no_ind_check_) {
      // inductive case check based on initial state predicates like in
      // Sheeran et al 2003: assert that s_0 is an initial state and no
      // other state s_1,...,s_{i} is an initial state + simple path
      // constraints.
      // NOTE: we do this check in a new push/pop frame, which does not
      // benefit from incrementality. At least, we should build a
      // conjunction of initial state constraints that is added back and
      // where we append a new conjunct for each time step.
      solver_->push();
      solver_->assert_formula(init0_);
      smt::Term neg_init_at_i = unroller_.at_time(
	solver_->make_term(Not, ts_.init()), i);
      neg_init_terms_ = solver_->make_term(And, neg_init_terms_, neg_init_at_i);
      solver_->assert_formula(neg_init_terms_);
      kind_log_msg(1, "", "checking inductive step (initial states) at bound: {}", i);
      res = solver_->check_sat_assuming(sel_assumption_);
      if (res.is_unsat()) {
	return ProverResult::TRUE;
      }
      solver_->pop();
    }

    // open new frame; this is to be able to remove bad state predicate added next
    solver_->push();

    // for inductive case and base case: add bad state predicate
    solver_->assert_formula(unroller_.at_time(bad_, i));

    // inductive case check
    if (!options_.kind_no_ind_check_) {
      kind_log_msg(1, "", "checking inductive step (property) at bound: {}", i);
      res = solver_->check_sat_assuming(sel_assumption_);
      if (res.is_unsat()) {
	return ProverResult::TRUE;
      }
    }

    // base case check

    //enable initial state predicate
    assert(sel_assumption_.size() == 1);
    sel_assumption_.pop_back();
    sel_assumption_.push_back(not_sel_init_);

    kind_log_msg(1, "", "checking base case at bound: {}", i);
    res = solver_->check_sat_assuming(sel_assumption_);
    if (res.is_sat()) {
      compute_witness();
      return ProverResult::FALSE;
    }

    solver_->pop();

    // add transition and negated bad state property
    // it is sound to add the negated bad state property for use in
    // next base case checks and inductive case checks (initial
    // states) because we proved in base check that it is implied when
    // assuming initial state predicate
    solver_->assert_formula(unroller_.at_time(ts_.trans(), i));
    solver_->assert_formula(unroller_.at_time(solver_->make_term(Not, bad_), i));

    reached_k_++;
  } //end: for all bounds
  
  return ProverResult::UNKNOWN;
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
  kind_log_msg(1, "", "checking simple path (eager) at bound: {}", i);

  const bool no_simp_path_check = options_.kind_no_simple_path_check_;

  //Here we assume that all pairs for values smaller than 'i' have been added
  //If simple path checking is disabled then we still need the final
  //solver call below for inductive case check
  for (int j = 0; (!no_simp_path_check && j < i); j++) {
    Term constraint = simple_path_constraint(j, i);
    kind_log_msg(3, "   ", "adding simple path clause for pair 'j,i' = {},{}", j,i);
    solver_->assert_formula(constraint);
  }

  kind_log_msg(2, "    ", "calling solver for simple path check");
  Result r = solver_->check_sat_assuming(sel_assumption_);
  if (r.is_unsat()) {
    kind_log_msg(2, "      ", "simple path check UNSAT");
    return true;
  }

  return false;
}

bool KInduction::check_simple_path_lazy(int i)
{
  kind_log_msg(1, "", "checking simple path (lazy) at bound: {}", i);
  bool added_to_simple_path = false;

  // If no_multi_call == true, then we add as many simple path
  // constraints as necessary based on the current model of the
  // solver. These constraints are collected in vector 'vec' first and
  // then added in one pass. Otherwise, the solver is called again
  // after *each* added constraint to get a new model.
  const bool no_multi_call = options_.kind_no_multi_call_simple_path_check_;
  smt::TermVec vec;
  
  do {
    assert(vec.size() == 0);
    kind_log_msg(2, "    ", "calling solver for simple path check");
    Result r = solver_->check_sat_assuming(sel_assumption_);
    if (r.is_unsat()) {
      kind_log_msg(2, "      ", "simple path check UNSAT");
      return true;
    }

    // We need the above check_sat call if simple path checks are
    // disabled because it is part of the inductive step check. Hence
    // we return immediately here and skip simple path checking.
    if (options_.kind_no_simple_path_check_)
      return false;

    added_to_simple_path = false;

    for (int j = 0; j < i && (no_multi_call || !added_to_simple_path); ++j) {
      for (int l = j + 1; l <= i; ++l) {
        Term constraint = simple_path_constraint(j, l);
	kind_log_msg(3, "    ", "checking constraint for pair j,l = {} , {}", j,l);
        if (solver_->get_value(constraint) == false_) {
	  kind_log_msg(3, "      ", "adding constraint for pair j,l = {} , {}", j,l);
          added_to_simple_path = true;
          if (!no_multi_call) {
            solver_->assert_formula(constraint);
            break;
          } else
            vec.push_back(constraint);
        }
      }
    }

    while (vec.size()) {
      assert(no_multi_call);
      Term constraint = vec.back();
      vec.pop_back();
      solver_->assert_formula(constraint);
    }

  } while (added_to_simple_path);

  return false;
}

template <typename... Args>
void KInduction::kind_log_msg(size_t level, const std::string & indent,
			      const std::string & format, const Args &... args)
{
  logger.log(level, indent + kind_engine_name_ + " " + format, args...);
}

}  // namespace pono
