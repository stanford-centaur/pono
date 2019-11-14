#include "bmc.h"
#include "utils/logger.h"

using namespace smt;

namespace cosa {

Bmc::Bmc(const Property & p, SmtSolver & solver) : super(p, solver)
{
  initialize();
}

Bmc::~Bmc() {}

void Bmc::initialize()
{
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
  for (int i = 0; i <= k; ++i) {
    if (!step(i)) {
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
    solver_->assert_formula(unroller_.at_time(ts_.trans(), i - 1));
  }

  solver_->push();
  logger.log(1, "Checking bmc at bound: {}", i);
  solver_->assert_formula(unroller_.at_time(bad_, i));
  Result r = solver_->check_sat();
  if (r.is_sat()) {
    res = false;
  } else {
    solver_->pop();
  }

  ++reached_k_;

  return res;
}

}  // namespace cosa
