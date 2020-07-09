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

#include "prover.h"
#include "available_solvers.h"

#include <climits>

using namespace smt;

namespace pono {

Prover::Prover(const Property & p, smt::SolverEnum se)
    : Prover(p, create_solver(se))
{
  solver_->set_opt("incremental", "true");
  solver_->set_opt("produce-models", "true");
}

Prover::Prover(const Property & p, smt::SmtSolver s)
    : solver_(s),
      to_prover_solver_(s),
      to_orig_ts_solver_(p.transition_system().solver()),
      property_(p, to_prover_solver_),
      ts_(property_.transition_system()),
      unroller_(ts_, solver_)
{
}

Prover::~Prover() {}

void Prover::initialize()
{
  reached_k_ = -1;
  bad_ = solver_->make_term(smt::PrimOp::Not, property_.prop());
}

ProverResult Prover::prove() { return check_until(INT_MAX); }

bool Prover::witness(std::vector<UnorderedTermMap> & out)
{
  // TODO: make sure the solver state is SAT

  for (int i = 0; i <= reached_k_; ++i) {
    out.push_back(UnorderedTermMap());
    UnorderedTermMap & map = out.back();

    for (auto v : ts_.statevars()) {
      Term vi = unroller_.at_time(v, i);
      Term r = solver_->get_value(vi);
      map[v] = r;
    }

    for (auto v : ts_.inputvars()) {
      Term vi = unroller_.at_time(v, i);
      Term r = solver_->get_value(vi);
      map[v] = r;
    }

    for (auto elem : ts_.named_terms()) {
      Term ti = unroller_.at_time(elem.second, i);
      map[elem.second] = solver_->get_value(ti);
    }
  }

  return true;
}

}  // namespace pono
