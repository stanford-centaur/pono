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
  bin_search_frames_ = 0;
  bound_step_ = opt.bmc_bound_step_;
  bound_start_ = opt.bmc_bound_start_;
}

Bmc::~Bmc() {}

void Bmc::initialize()
{
  if (initialized_) {
    return;
  }

  super::initialize();

  // NOTE: for any engine; There's an implicit assumption that this solver is only used for
  // model checking once Otherwise there could be conflicting assertions to
  // the solver or it could just be polluted with redundant assertions in the
  // future we can use solver_->reset_assertions(), but it is not currently
  // supported in boolector
  logger.log(2, "DEBUG adding init constraint for step 0");
  solver_->assert_formula(unroller_.at_time(ts_.init(), 0));
}

ProverResult Bmc::check_until(int k)
{
  initialize();

  //NOTE/TODO: there is a corner case where an instance is trivially
  //unsatisfiable, i.e., safe, when the conjunction of initial state
  //predicate and transition (+ any constraints) is already unsat. We
  //could also check this using unsat core functionality of solver (if
  //supported), and check if bad state predicate is in core

  logger.log(1, "DEBUG BMC bound_start_ {} ", bound_start_);  
  logger.log(1, "DEBUG BMC bound_step_ {} ", bound_step_);  
  
  // reached_k == -1 initially
  
  //  if (start_bound > 0)
  //  reached_k_ = start_bound - 1;
  
//  for (int j = 1; j < start_bound; ++j)
  // {
  //   std::cout << "DEBUG (check-until) adding trans for j-1 == " << j - 1 << std::endl;
  //   solver_->assert_formula(unroller_.at_time(ts_.trans(), j - 1));
//  }

  const bool exp_step = options_.bmc_exponential_step_;

  for (int i = bound_start_; i <= k; i = exp_step ? (i == 0 ? 1 : i << 1) : (i + bound_step_)) {
    if (!step(i)) {
      compute_witness();
      return ProverResult::FALSE;
    }
  }
  return ProverResult::UNKNOWN;
}

bool Bmc::step(int i)
{
  logger.log(1, "Checking bmc at bound: {}", i);
  
  if (i <= reached_k_) {
    return true;
  }

  bool res = true;
//CHECK code simplification: if-statement needed? seems not needed
  if (i > 0) {
//
    logger.log(2, "DEBUG reached k {}, i {} ", reached_k_, i);
    for (int j = reached_k_ == -1 ? 1 : reached_k_ + 1; j <= i; j++) {
      logger.log(2, "DEBUG adding trans for j-1 == {}", j - 1);
      solver_->assert_formula(unroller_.at_time(ts_.trans(), j - 1));
      if (options_.bmc_neg_init_step_) {
	logger.log(2, "DEBUG adding negated init constraint for step {}", j);
	Term not_init = solver_->make_term(PrimOp::Not, unroller_.at_time(ts_.init(), j));
	solver_->assert_formula(not_init);
      }
    }
    //OLD
//    solver_->assert_formula(unroller_.at_time(ts_.trans(), i - 1));
  }

  solver_->push();

  const int cex_guarantee = !options_.bmc_single_bad_state_;

  Term clause;
  if (cex_guarantee) {
    // make sure we cover *all* states
  // TODO (not critical): can make 'solver_->make_term(false)' a constant in the object
    clause = solver_->make_term(false);
    for (int j = reached_k_ + 1; j <= i; j++) {
      logger.log(2, "DEBUG adding bad state constraint for j == {}", j);
      clause = solver_->make_term(PrimOp::Or, clause, unroller_.at_time(bad_, j));
    }
  } else {
    logger.log(2, "DEBUG adding bad state constraint for i == {}", i);
    clause = unroller_.at_time(bad_, i);
  }
  
  solver_->assert_formula(clause);
  
  //TODO: add bad state clause here
    // OLD
    // solver_->assert_formula(unroller_.at_time(bad_, i));


  Result r = solver_->check_sat();
  if (r.is_sat()) {
    logger.log(1, "  bmc check at bound {} satisfiable", i);
    res = false;
    if (cex_guarantee) {
      logger.log(2, "DEBUG saving reached_k_ = {}", reached_k_);
      int reached_k_saved = reached_k_;
      int cex_upper_bound = bmc_interval_get_cex_ub(reached_k_ + 1, i);
      //FIX for corner case of length-1 interval: some solvers don't
      //allow 'push' before 'get-value' even when no terms are asserted
      //before 'get-value'. Hence don't 'push' if we don't add any terms
      //in function 'bmc_interval_block_cex_ub'.
      if (cex_upper_bound + 1 <= i)
	solver_->push();
      bmc_interval_block_cex_ub(cex_upper_bound + 1, i);
      // find shortest cex within tested interval given by bad state clause
//    bmc_interval_find_shortest_cex(cex_upper_bound);
      bool success = bmc_interval_find_shortest_cex_binary_search(cex_upper_bound);
      if (!success) {
	reached_k_ = reached_k_saved;
	//clear constraints added during upper bound computation and binary search
	if (cex_upper_bound + 1 <= i)
	  solver_->pop();
	while(bin_search_frames_-- > 0)
	  solver_->pop();
	bmc_interval_find_shortest_cex_linear_search(cex_upper_bound);
      }
    } else {
      // Handle corner case when using single bad state constraints and
      // interval search: for witness printing, which depends on
      // reached_k_, we must set reached_k_ to the bound that preceeds
      // the bound 'i' where the cex was found
      reached_k_ = i - 1;
    }
  } else {
    logger.log(1, "  bmc check at bound {} unsatisfiable", i);
    solver_->pop();
    // TODO: could add expert option to add *all* since last reached_k,
    // whether being tested above or not
    if (options_.bmc_neg_bad_step_) {
      Term not_bad;
      if (cex_guarantee) {
	for (int j = reached_k_ + 1; j <= i; j++) {
	  logger.log(2, "DEBUG adding negated bad state constraint for j == {}", j);
	  not_bad = solver_->make_term(PrimOp::Not, unroller_.at_time(bad_, j));
	  solver_->assert_formula(not_bad);
	}
      } else {
	logger.log(2, "DEBUG adding negated bad state constraint for i == {}", i);
	not_bad = solver_->make_term(PrimOp::Not, unroller_.at_time(bad_, i));
	solver_->assert_formula(not_bad);
      }
    }
    reached_k_ = i;
  }

  return res;
}
  
int Bmc::bmc_interval_get_cex_ub(const int lb, const int ub)
{
  const Term true_term = solver_->make_term(true);
  assert(lb <= ub);
  
  logger.log(2, "DEBUG get cex upper bound: lower bound = {}, upper bound = {} ", lb, ub);
  
  int j;
  for (j = lb; j <= ub; j++) {
    Term bad_state_at_j = unroller_.at_time(bad_, j);
    logger.log(2, "DEBUG get cex upper bound, checking value of bad state constraint j = {}", j);
    //TODO check: proper use of return value of "get_value"?
    if (solver_->get_value(bad_state_at_j) == true_term) {
      logger.log(2, "DEBUG get cex upper bound, found at j = {}", j);
      break;
    }
  }
  assert(j <= ub);
  
  return j;
}

int Bmc::bmc_interval_block_cex_ub(const int start, const int end)
{
  //TODO CHECK that witness printing still works, i.e., we don't add new constraints after a sat-call
  //for search of shortest cex: block bad state constraints in interval [start,end]
  logger.log(2, "DEBUG get cex upper bound, permanently blocking [start,end] = [{},{}]", start, end); 
  for (int k = start; k <= end; k++) {
    Term not_bad = solver_->make_term(PrimOp::Not, unroller_.at_time(bad_, k));
    logger.log(3, "DEBUG get cex upper bound, "	\
	       "adding permanent blocking bad state constraint for k == {}", k);
    solver_->assert_formula(not_bad);
  }
}
  
bool Bmc::bmc_interval_find_shortest_cex_binary_search(const int upper_bound)
{
  assert (bin_search_frames_ == 0);
  assert (reached_k_ < upper_bound);
  logger.log(2, "DEBUG binary search, cex found in interval "\
	     "[reached_k+1,upper_bound] = [{},{}]", reached_k_ + 1, upper_bound);

  if (upper_bound - reached_k_ == 1) {
    logger.log(2, "DEBUG interval has length 1, skipping search for shortest cex");
    return true;
  }

  int low = reached_k_ + 1;
  int high = upper_bound;
  while (low <= high) {
    //Term clause = solver_->make_term(false);
    int mid = low + (high - low) / 2;
    logger.log(2, "\nDEBUG binary search, (low, mid, high) = ({}, {}, {})", low, mid, high);

    logger.log(3, "DEBUG binary search, solver->push()");
    //solver_->pop();
    solver_->push();
    bin_search_frames_++;

//TODO REMOVE
    //  assert(solver_->check_sat().is_sat());

    
    int j;
//    for (j = low; j <= mid; j++) {
    //we search for cex in [low,mid] hence block [mid+1,high]
    logger.log(2, "DEBUG binary search, searching for cex in [low,mid] = [{},{}]", low, mid);
    logger.log(2, "DEBUG binary search, temporarily blocking [mid+1,high] = [{},{}]", mid + 1, high); 
    for (j = mid + 1; j <= high; j++) {
      logger.log(3, "DEBUG binary search, finding shortest cex---"\
		 "adding blocking bad state constraint for j == {}", j);
      //clause = solver_->make_term(PrimOp::Or, clause, unroller_.at_time(bad_, j));
      Term not_bad = solver_->make_term(PrimOp::Not, unroller_.at_time(bad_, j));
      solver_->assert_formula(not_bad);
    }
//    solver_->assert_formula(clause);

    Result r = solver_->check_sat();
    assert(r.is_sat() || r.is_unsat());
    if (r.is_sat()) {
      logger.log(2, "DEBUG binary search, sat result: {}", r);
      logger.log(2, "DEBUG binary search, cex found in [low,mid] = [{},{}]", low, mid);
      logger.log(2, "DEBUG binary search, [mid+1,high] = [{},{}] now permanently blocked", mid + 1, high); 
      // if low == mid in current iteration, then we have tested a single
      // bad state constraint; can exit loop in case of satisfiability
      if (low == mid)
	break;
      else {
	high = bmc_interval_get_cex_ub(low, mid);
	bmc_interval_block_cex_ub(high + 1, mid);
      }
    } else {
      logger.log(2, "DEBUG binary search, unsat result: {}", r);
      logger.log(2, "DEBUG binary search, no cex in [low,mid] = [{},{}]", low, mid);
      if (low >= high) {
	// handle rare corner case
	logger.log(1, "BMC FAILURE: formula overconstrained,"\
		      " falling back to linear search");
	return false;
      }
      logger.log(2, "DEBUG binary search, unblocking [mid+1,high] = [{},{}]", mid + 1, high);
      //remove previoulsy added blocking constraints for [mid+1,high]
      logger.log(3, "DEBUG binary search, solver->pop()");
      solver_->pop();
      assert(bin_search_frames_ > 0);
      bin_search_frames_--;
      // update reached k; we have iteratively shown that no cex exists in [0,mid]
      reached_k_ = mid;
      logger.log(2, "DEBUG binary search, permanently blocking [low,mid] = [{},{}]", low, mid); 
      //no cex found in [low,mid] hence block [low,mid]
      for (j = low; j <= mid; j++) {
	logger.log(3, "DEBUG binary search, finding shortest cex---"	\
		   "adding blocking bad state constraint for j == {}", j);
	Term not_bad = solver_->make_term(PrimOp::Not, unroller_.at_time(bad_, j));
	solver_->assert_formula(not_bad);
      }

//TODO REMOVE
      //  assert(solver_->check_sat().is_sat());
      
      low = mid + 1;
    }
  }

  //must find cex inside sat-branch of if-then-else above
  assert(low <= high);
  //reached_k_ has been correctly updated to low - 1, i.e., cex bound - 1
  assert(reached_k_ + 1 == low);
  logger.log(1, "DEBUG binary search, shortest cex at bound low == {},"\
	     " reached_k = {}", low, reached_k_);
  return true;
}
  
void Bmc::bmc_interval_find_shortest_cex_linear_search(const int upper_bound)
{  
  assert (reached_k_ < upper_bound);
  logger.log(2, "DEBUG linear search, cex found in interval: lower bound = reached k = {},"\
	     " upper bound = {}", reached_k_, upper_bound);

//TODO: immediately return for length-1 intervals

  if (upper_bound - reached_k_ == 1) {
    logger.log(2, "DEBUG interval has length 1, skipping search for shortest cex");
    return;
  }
  

//TODO: below loop searches iteratively; could do exponential
//steps also, like in "check_until"; BMC with exponential steps could
//be useful if there is a very deep bug at a large bound and the costs
//of the SAT calls does not increase very much with the increased
//bound, but still the cumulative time for the SAT calls at all bounds
//is higher
  
//  return;
  int j;
  for (j = reached_k_ + 1; j <= upper_bound; j++) {
    //pop: remove the latest bad state clause
    solver_->pop();
    solver_->push();
    logger.log(2, "DEBUG finding shortest cex---adding bad state constraint for j == {}", j);
    solver_->assert_formula(unroller_.at_time(bad_, j));
    if (solver_->check_sat().is_sat()) {
      break;
    }
    else
      reached_k_ = j;
  }
  // must have found cex in the interval
  if (j > upper_bound)
    throw PonoException("BMC FAILURE in linear search: formula overconstrained");

  assert(reached_k_ + 1 == j);
  logger.log(1, "DEBUG linear search found shortest cex at bound j == {}, reached_k {}", j, reached_k_);
}
  
}  // namespace pono
