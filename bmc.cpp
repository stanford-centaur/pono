#include "bmc.h"

using namespace smt;

namespace cosa
{

  Bmc::Bmc(const Property &p, SmtSolver &solver):
    ts_(p.transition_system()),
    property_(p),
    solver_(solver),
    unroller_(ts_, solver_)
  {
    initialize();
  }

  Bmc::~Bmc()
  {
  }

  void Bmc::initialize()
  {
    reached_k_ = -1;
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

  ProverResult Bmc::prove()
  {
    for (int i = 0; ; ++i) {
      if (!step(i)) {
	return ProverResult::FALSE;
      }
    }
    return ProverResult::UNKNOWN;
  }

  bool Bmc::witness(std::vector<UnorderedTermMap> &out)
  {
    // TODO: make sure the solver state is SAT

    for (int i = 0; i <= reached_k_; ++i) {
      out.push_back(UnorderedTermMap());
      UnorderedTermMap &map = out.back();

      for (auto v : ts_.states()) {
	Term vi = unroller_.at_time(v, i);
	Term r = solver_->get_value(vi);
	map[v] = r;
      }

      if (i != reached_k_) {
	for (auto v : ts_.inputs()) {
	  Term vi = unroller_.at_time(v, i);
	  Term r = solver_->get_value(vi);
	  map[v] = r;
	}
      }
    }

    return true;
  }

  bool Bmc::step(int i)
  {
    if (i <= reached_k_) {
      return true;
    }

    std::cout << "Checking BMC Bound " << i << std::endl;

    bool res = true;
    Term bad = solver_->make_term(PrimOp::Not, property_.prop());

    if (i > 0) {
      solver_->assert_formula(unroller_.at_time(ts_.trans(), i-1));
    }

    solver_->push();
    solver_->assert_formula(unroller_.at_time(bad, i));
    Result r = solver_->check_sat();
    if (r.is_sat()) {
      res = false;
    } else {
      solver_->pop();
    }

    ++reached_k_;

    return res;
  }

} // namespace cosa
