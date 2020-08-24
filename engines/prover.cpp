/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann, Florian Lonsing
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
#include <functional>

using namespace smt;
using namespace std;

namespace pono {

Prover::Prover(const Property & p, smt::SolverEnum se)
    : Prover(p, create_solver(se))
{
  solver_->set_opt("incremental", "true");
  solver_->set_opt("produce-models", "true");
}

Prover::Prover(const Property & p, const smt::SmtSolver & s)
    : solver_(s),
      to_prover_solver_(s),
      property_(p, to_prover_solver_),
      ts_(property_.transition_system()),
      orig_ts_(p.transition_system()),
      unroller_(ts_, solver_)
{
}

Prover::Prover(const PonoOptions & opt, const Property & p, smt::SolverEnum se)
    : Prover(opt, p, create_solver(se))
{
  solver_->set_opt("incremental", "true");
  solver_->set_opt("produce-models", "true");
}

Prover::Prover(const PonoOptions & opt,
               const Property & p,
               const smt::SmtSolver & s)
    : solver_(s),
      to_prover_solver_(solver_),
      property_(p, to_prover_solver_),
      ts_(property_.transition_system()),
      orig_ts_(p.transition_system()),
      unroller_(ts_, solver_),
      options_(opt)
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

  function<Term(const Term &, SortKind)> transfer_to_prover_as;
  function<Term(const Term &, SortKind)> transfer_to_orig_ts_as;
  TermTranslator to_orig_ts_solver(orig_ts_.solver());
  if (solver_ == orig_ts_.solver()) {
    // don't need to transfer terms if the solvers are the same
    transfer_to_prover_as = [](const Term & t, SortKind sk) { return t; };
    transfer_to_orig_ts_as = [](const Term & t, SortKind sk) { return t; };
  } else {
    // need to add symbols to cache
    UnorderedTermMap & cache = to_orig_ts_solver.get_cache();
    for (auto v : orig_ts_.statevars()) {
      cache[to_prover_solver_.transfer_term(v)] = v;
    }
    for (auto v : orig_ts_.inputvars()) {
      cache[to_prover_solver_.transfer_term(v)] = v;
    }

    transfer_to_prover_as = [this](const Term & t, SortKind sk) {
      return to_prover_solver_.transfer_term(t, sk);
    };
    transfer_to_orig_ts_as = [&to_orig_ts_solver](const Term & t, SortKind sk) {
      return to_orig_ts_solver.transfer_term(t, sk);
    };
  }

  for (int i = 0; i <= reached_k_; ++i) {
    out.push_back(UnorderedTermMap());
    UnorderedTermMap & map = out.back();

    for (auto v : orig_ts_.statevars()) {
      SortKind sk = v->get_sort()->get_sort_kind();
      Term vi = unroller_.at_time(transfer_to_prover_as(v, sk), i);
      Term r = solver_->get_value(vi);
      map[v] = transfer_to_orig_ts_as(r, sk);
    }

    for (auto v : orig_ts_.inputvars()) {
      SortKind sk = v->get_sort()->get_sort_kind();
      Term vi = unroller_.at_time(transfer_to_prover_as(v, sk), i);
      Term r = solver_->get_value(vi);
      map[v] = transfer_to_orig_ts_as(r, sk);
    }

    for (auto elem : orig_ts_.named_terms()) {
      SortKind sk = elem.second->get_sort()->get_sort_kind();
      Term ti = unroller_.at_time(transfer_to_prover_as(elem.second, sk), i);
      map[elem.second] = transfer_to_orig_ts_as(solver_->get_value(ti), sk);
    }
  }

  return true;
}

}  // namespace pono
