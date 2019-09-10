#include "kinduction.h"

using namespace smt;

namespace cosa
{

  KInduction::KInduction(const Property &p, smt::SmtSolver &solver):
    ts_(p.transition_system()),
    property_(p),
    solver_(solver),
    unroller_(ts_, solver_)
  {
    initialize();
  }

  KInduction::~KInduction()
  {
  }

  KInduction::initialize()
  {
    reached_k_ = -1;
    //solver_->reset_assertions();
  }

  bool KInduction::check_until(size_t k)
  {
    for (size_t i = 0; i <= k; ++i) {
      if (!base_step(i)) {
	return false;
      }
      if (inductive_step(i)) {
	return true;
      }
    }
    //TODO: return unknown
  }
  
  bool KInduction::base_step(size_t i)
  {
    if (i <= reached_k_) {
      return true;
    }

    Term bad = solver_->make_term(PrimOp::Not, property_.prop());
    
    solver_->push();
    solver_->assert_formula(unroller_.at_time(ts_.init(), 0));
    solver_->assert_formula(unroller_.at_time(bad, i));
    Result r = solver_->check_sat();
    if (r.is_sat) {
      return false;
    }
    solver_->pop();

    solver_->assert_formula(unroller_.at_time(ts_.trans(), i));
    solver_->assert_formula(unroller_.at_time(property_.prop, i));
    
    return true;
  }

  bool KInduction::inductive_step(size_t i)
  {
    if (i <= reached_k_) {
      return false;
    }

    Term bad = solver_->make_term(PrimOp::Not, property_.prop());

    solver_->push();
    solver_->assert_formula(unroller_.at_time(bad, i+1));
    Result r = solver_->check_sat();
    if (r.is_unsat) {
      return true;
    }
    solver_->pop();

    ++reached_k_;

    return false;
  }

  void KInduction::add_simple_path_constraint(size_t i, size_t j)
  {
    //TODO
  }
  
} // namespace cosa
