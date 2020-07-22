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
  // TODO: Figure out if we need this
  //       shouldn't have to do this special check
  if (reached_k_ < 1) {
    solver_->push();
    solver_->assert_formula(ts_.init());
    solver_->assert_formula(bad_);
    Result r = solver_->check_sat();
    if (r.is_sat()) {
      return ProverResult::FALSE;
    } else {
      assert(r.is_unsat());
      reached_k_ = 1;  // keep reached_k_ aligned with number of frames
    }
    solver_->pop();
  }

  while (reached_k_ <= k) {
    assert(reached_k_ == frames_.size());
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

bool ModelBasedIC3::intersects_bad()
{
  solver_->push();

  // assert the frame (conjunction over clauses)
  for (auto c : frames_.back()) {
    solver_->assert_formula(c);
  }

  // see if it intersects with bad
  solver_->assert_formula(bad_);

  Result r = solver_->check_sat();
  if (r.is_sat()) {
    // create a proof goal for the bad state
    const UnorderedTermSet & statevars = ts_.statevars();
    TermVec cube_vec;
    cube_vec.reserve(statevars.size());
    Term eq;
    for (auto sv : statevars) {
      eq = solver_->make_term(Equal, sv, solver_->get_value(sv));
      cube_vec.push_back(eq);
    }
    Cube c(solver_, cube_vec);
    proof_goals_[reached_k_].push_back(c);
  }

  solver_->pop();

  assert(!r.is_unknown());
  return r.is_sat();
}

}  // namespace pono
