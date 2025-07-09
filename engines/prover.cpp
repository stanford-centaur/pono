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
#include <cstddef>
#include <exception>
#include <functional>
#include <vector>

#include "core/ts.h"
#include "options/options.h"
#include "smt-switch/smt.h"

using namespace smt;
using namespace std;

namespace pono {
BaseProver::BaseProver(const TransitionSystem & ts, const smt::SmtSolver & s)
    : BaseProver(ts, s, PonoOptions())
{
}

BaseProver::BaseProver(const TransitionSystem & ts,
                       const smt::SmtSolver & s,
                       PonoOptions opt)
    : solver_(s),
      to_prover_solver_(s),
      orig_ts_(ts),
      ts_(ts, to_prover_solver_),
      options_(opt)
{
}

void BaseProver::initialize() { initialized_ = true; }

ProverResult BaseProver::prove() { return check_until(INT_MAX); }

Term BaseProver::to_orig_ts(Term t, SortKind sk)
{
  if (solver_ == orig_ts_.solver()) {
    // don't need to transfer terms if the solvers are the same
    return t;
  } else {
    // need to add symbols to cache
    TermTranslator to_orig_ts_solver(orig_ts_.solver());
    UnorderedTermMap & cache = to_orig_ts_solver.get_cache();
    for (const auto & v : orig_ts_.statevars()) {
      cache[to_prover_solver_.transfer_term(v)] = v;
      const Term & nv = orig_ts_.next(v);
      cache[to_prover_solver_.transfer_term(nv)] = v;
    }
    for (const auto & v : orig_ts_.inputvars()) {
      cache[to_prover_solver_.transfer_term(v)] = v;
    }
    // TODO: need a to add UFs to the cache also
    return to_orig_ts_solver.transfer_term(t, sk);
  }
}

Term BaseProver::to_orig_ts(Term t)
{
  return to_orig_ts(t, t->get_sort()->get_sort_kind());
}

TransitionSystem & BaseProver::prover_interface_ts() { return ts_; }

SafetyProver::SafetyProver(const SafetyProperty & p,
                           const TransitionSystem & ts,
                           const smt::SmtSolver & s,
                           PonoOptions opt)
    : BaseProver(ts, s, opt),
      orig_property_(p),
      unroller_(ts_),
      bad_(solver_->make_term(
          smt::PrimOp::Not,
          ts_.solver() == orig_property_.solver()
              ? orig_property_.prop()
              : to_prover_solver_.transfer_term(orig_property_.prop(), BOOL)))
{
}

void SafetyProver::initialize()
{
  if (initialized_) {
    return;
  }

  if (!ts_.only_curr(bad_)) {
    throw PonoException(
        "Safety property should not contain inputs or next state variables");
  }

  BaseProver::initialize();
}

bool SafetyProver::witness(std::vector<UnorderedTermMap> & out)
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
    // need to add symbols to cache
    UnorderedTermMap & cache = to_orig_ts_solver.get_cache();
    for (const auto & v : orig_ts_.statevars()) {
      cache[to_prover_solver_.transfer_term(v)] = v;
    }
    for (const auto & v : orig_ts_.inputvars()) {
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

    for (const auto & v : orig_ts_.statevars()) {
      const SortKind & sk = v->get_sort()->get_sort_kind();
      const Term & pv = transfer_to_prover_as(v, sk);
      map[v] = transfer_to_orig_ts_as(wit_map.at(pv), sk);
    }

    for (const auto & v : orig_ts_.inputvars()) {
      const SortKind & sk = v->get_sort()->get_sort_kind();
      const Term & pv = transfer_to_prover_as(v, sk);
      try {
        map[v] = transfer_to_orig_ts_as(wit_map.at(pv), sk);
      }
      catch (std::exception & e) {
        success = false;
        break;
      }
    }

    if (success) {
      for (const auto & elem : orig_ts_.named_terms()) {
        const SortKind & sk = elem.second->get_sort()->get_sort_kind();
        const Term & pt = transfer_to_prover_as(elem.second, sk);
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

size_t SafetyProver::witness_length() const { return reached_k_ + 1; }

Term SafetyProver::invar()
{
  if (!invar_) {
    throw PonoException(
        "Failed to return invar. Be sure that the property was proven by an "
        "engine the supports returning invariants.");
  }
  return to_orig_ts(invar_, BOOL);
}

bool SafetyProver::compute_witness()
{
  // TODO: make sure the solver state is SAT

  for (size_t i = 0, wl = witness_length(); i <= wl; ++i) {
    witness_.push_back(UnorderedTermMap());
    UnorderedTermMap & map = witness_.back();

    for (const auto & v : ts_.statevars()) {
      const Term & vi = unroller_.at_time(v, i);
      const Term & r = solver_->get_value(vi);
      map[v] = r;
    }

    for (const auto & v : ts_.inputvars()) {
      const Term & vi = unroller_.at_time(v, i);
      const Term & r = solver_->get_value(vi);
      map[v] = r;
    }

    for (const auto & elem : ts_.named_terms()) {
      const Term & ti = unroller_.at_time(elem.second, i);
      map[elem.second] = solver_->get_value(ti);
    }
  }

  return true;
}

}  // namespace pono
