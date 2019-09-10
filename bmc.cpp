#include "bmc.h"

using namespace smt;

namespace cosa
{

  Bmc::Bmc(const RelationalTransitionSystem &ts, const Property &p,
	   SmtSolver &solver):
    ts_(ts),
    property_(p),
    solver_(solver),
    unroller_(ts, solver)
  {
    solver_->assert_formula(unroller_.at_time(ts_.init(), 0));
  }

  Bmc::~Bmc()
  {
  }

  bool Bmc::check_until(size_t k)
  {
    size_t i = 0;
    while (i <= k) {
      if (!step(i)) {
	return false;
      }
    }
    return true;
  }

  bool Bmc::step(size_t i)
  {
    bool res = true;
    Term bad = solver_->make_term(PrimOp::Not, property_.prop());

    if (i > 0) {
      solver_->assert_formula(unroller_.at_time(ts_.trans(), i-1));
    }

    //solver_->push();
    solver_->assert_formula(unroller_.at_time(bad, i));
    Result r = solver_->check_sat();
    if (r.is_sat()) {
      res = false;
    } else {
      //solver_->pop();
    }
    
    return res;
  }
  
} // namespace cosa
