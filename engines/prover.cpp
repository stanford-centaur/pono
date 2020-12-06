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

#include "engines/prover.h"

#include <cassert>
#include <climits>
#include <functional>

#include "modifiers/coi.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

Prover::Prover(Property & p, smt::SolverEnum se)
    : Prover(p, create_solver(se))
{
  solver_->set_opt("incremental", "true");
  solver_->set_opt("produce-models", "true");
}

Prover::Prover(Property & p, const smt::SmtSolver & s)
    : solver_(s),
      to_prover_solver_(s),
      property_(p, to_prover_solver_),
      ts_(&property_.transition_system()),
      orig_ts_(p.transition_system()),
      unroller_(*ts_, solver_)
{
}

Prover::Prover(const PonoOptions & opt, Property & p, smt::SolverEnum se)
    : Prover(opt, p, create_solver(se))
{
  solver_->set_opt("incremental", "true");
  solver_->set_opt("produce-models", "true");
}

Prover::Prover(const PonoOptions & opt, Property & p, const smt::SmtSolver & s)
    : solver_(s),
      to_prover_solver_(solver_),
      property_(p, to_prover_solver_),
      ts_(&property_.transition_system()),
      orig_ts_(p.transition_system()),
      unroller_(*ts_, solver_),
      options_(opt)
{
}

Prover::~Prover() {}

void Prover::initialize()
{
  reached_k_ = -1;
  bad_ = solver_->make_term(smt::PrimOp::Not, property_.prop());
  assert(ts_->only_curr(bad_));
  if (options_.static_coi_) {
    /* Compute the set of state/input variables related to the
       bad-state property. Based on that information, rebuild the
       transition relation of the transition system. */
    ConeOfInfluence coi(*ts_, { bad_ }, {}, options_.verbosity_);
  }
}

ProverResult Prover::prove() { return check_until(INT_MAX); }

bool Prover::witness(std::vector<UnorderedTermMap> & out)
{
  if (!witness_.size()) {
    throw PonoException(
        "Recovering witness failed. Make sure that there was "
        "a counterexample and that the engine supports witness generation.");
  }

  function<Term(const Term &, SortKind)> transfer_to_prover_as;
  function<Term(const Term &, SortKind)> transfer_to_orig_ts_as;
  TermTranslator to_orig_ts_solver(orig_ts_.solver());
  if (solver_ == orig_ts_.solver()) {
    // don't need to transfer terms if the solvers are the same
    transfer_to_prover_as = [](const Term & t, SortKind sk) { return t; };
    transfer_to_orig_ts_as = [](const Term & t, SortKind sk) { return t; };
  } else {
    /* TODO: double-check that transferring terms still works as
       intended in this branch when COI is used. */
    if (options_.static_coi_)
      throw PonoException(
          "Temporary restriction: cone-of-influence analysis "
          "currently incompatible with witness generation.");
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

  bool success = true;

  // Some backends don't support full witnesses
  // it will still populate state variables, but will return false instead of
  // true
  for (auto wit_map : witness_) {
    out.push_back(UnorderedTermMap());
    UnorderedTermMap & map = out.back();

    for (auto v : orig_ts_.statevars()) {
      SortKind sk = v->get_sort()->get_sort_kind();
      Term pv = transfer_to_prover_as(v, sk);
      map[v] = transfer_to_orig_ts_as(wit_map.at(pv), sk);
    }

    for (auto v : orig_ts_.inputvars()) {
      SortKind sk = v->get_sort()->get_sort_kind();
      Term pv = transfer_to_prover_as(v, sk);
      try {
        map[v] = transfer_to_orig_ts_as(wit_map.at(pv), sk);
      }
      catch (std::exception & e) {
        success = false;
        break;
      }
    }

    if (success) {
      for (auto elem : orig_ts_.named_terms()) {
        SortKind sk = elem.second->get_sort()->get_sort_kind();
        Term pt = transfer_to_prover_as(elem.second, sk);
        try {
          map[elem.second] = transfer_to_orig_ts_as(wit_map.at(pt), sk);
        }
        catch (std::exception & e) {
          success = false;
          break;
        }
      }
    }
  }

  return success;
}

Term Prover::invar()
{
  if (!invar_)
  {
    throw PonoException("Failed to return invar. Be sure that the property was proven "
                        "by an engine the supports returning invariants.");
  }
  return to_orig_ts(invar_, BOOL);
}

Term Prover::to_orig_ts(Term t, SortKind sk)
{
  if (solver_ == orig_ts_.solver()) {
    // don't need to transfer terms if the solvers are the same
    return t;
  } else {
    /* TODO: double-check that transferring terms still works as
       intended in this branch when COI is used. */
    if (options_.static_coi_)
      throw PonoException(
          "Temporary restriction: cone-of-influence analysis "
          "currently incompatible with witness generation.");
    // need to add symbols to cache
    TermTranslator to_orig_ts_solver(orig_ts_.solver());
    UnorderedTermMap & cache = to_orig_ts_solver.get_cache();
    for (auto v : orig_ts_.statevars()) {
      cache[to_prover_solver_.transfer_term(v)] = v;
      Term nv = orig_ts_.next(v);
      cache[to_prover_solver_.transfer_term(nv)] = v;
    }
    for (auto v : orig_ts_.inputvars()) {
      cache[to_prover_solver_.transfer_term(v)] = v;
    }
    return to_orig_ts_solver.transfer_term(t, sk);
  }
}

Term Prover::to_orig_ts(Term t)
{
  return to_orig_ts(t, t->get_sort()->get_sort_kind());
}

bool Prover::compute_witness()
{
  // TODO: make sure the solver state is SAT

  for (int i = 0; i <= reached_k_; ++i) {
    witness_.push_back(UnorderedTermMap());
    UnorderedTermMap & map = witness_.back();

    for (auto v : ts_->statevars()) {
      Term vi = unroller_.at_time(v, i);
      Term r = solver_->get_value(vi);
      map[v] = r;
    }

    for (auto v : ts_->inputvars()) {
      Term vi = unroller_.at_time(v, i);
      Term r = solver_->get_value(vi);
      map[v] = r;
    }

    for (auto elem : ts_->named_terms()) {
      Term ti = unroller_.at_time(elem.second, i);
      map[elem.second] = solver_->get_value(ti);
    }
  }

  return true;
}

}  // namespace pono
