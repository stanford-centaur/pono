#include "mbic3.h"

namespace pono {

ModelBasedIC3::ModelBasedIC3(const Property & p, const smt::SmtSolver & slv)
    : super(p, slv),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false)),
{
  initialize();
}

ModelBasedIC3::ModelBasedIC3(const PonoOptions & opt,
                             const Property p,
                             smt::SmtSolver slv)
    : super(opt, p, slv),
      true_(solver_->make_term(true)),
      false_(solver_->make_term(false)),
{
  initialize();
}

void ModelBasedIC3::initialize()
{
  super::initialize();
  // first frame is always the initial states
  frames_ = { ts_.init() };
  push_frame();
}

ProverResult ModelBasedIC3::check_until(int k)
{
  while (reached_k_ <= k) {
    // blocking phase
    while (intersects_bad()) {
      assert(!proof_goals_.empty());
      if (!block_all()) {
        // counter-example
        return ProverResult::FALSE;
      }
    }

    // propagation phase
    push_frame();
    propagate();
    if (is_proven()) {
      return ProverResult::TRUE;
    }
  }

  return ProverResult::UNKNOWN;
}

}  // namespace pono
