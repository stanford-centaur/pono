/*********************                                                        */
/*! \file msat_ic3ia.h
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief A backend using the open-source IC3IA implementation from
**        Alberto Griggio. Available here:
**        https://es-static.fbk.eu/people/griggio/ic3ia/index.html
**
**/

#include "engines/msat_ic3ia.h"

#include <cassert>

#include "smt-switch/msat_solver.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

MsatIC3IA::MsatIC3IA(const SafetyProperty & p,
                     const TransitionSystem & ts,
                     const SmtSolver & solver,
                     PonoOptions opt)
    : super(p, ts, solver, opt)
{
  engine_ = Engine::MSAT_IC3IA;
}

ProverResult MsatIC3IA::prove()
{
  initialize();

  if (ts_.solver()->get_solver_enum() != MSAT) {
    throw PonoException("MsatIC3IA only supports mathsat solver.");
  }

  // move everything over to a fresh solver
  shared_ptr<MsatSolver> msat_solver = static_pointer_cast<MsatSolver>(
      create_solver_for(MSAT, MSAT_IC3IA, false));
  msat_env env = msat_solver->get_msat_env();
  ::ic3ia::TransitionSystem ic3ia_ts(env);

  TermTranslator to_msat_solver(msat_solver);
  TermTranslator to_ts_solver(solver_);

  // give mapping between symbols
  UnorderedTermMap & ts_solver_cache = to_ts_solver.get_cache();
  for (const auto & v : ts_.statevars()) {
    ts_solver_cache[to_msat_solver.transfer_term(v)] = v;
    ts_solver_cache[to_msat_solver.transfer_term(ts_.next(v))] = ts_.next(v);
  }
  for (const auto & v : ts_.inputvars()) {
    ts_solver_cache[to_msat_solver.transfer_term(v)] = v;
  }

  // need to handle UFs also
  UnorderedTermSet ufs;
  auto is_uf = [](const Term & term) {
    return term->get_sort()->get_sort_kind() == smt::FUNCTION;
  };
  get_matching_terms(ts_.init(), ufs, is_uf);
  get_matching_terms(ts_.trans(), ufs, is_uf);
  get_matching_terms(bad_, ufs, is_uf);

  for (const auto & uf : ufs) {
    assert(uf->get_sort()->get_sort_kind() == smt::FUNCTION);
    ts_solver_cache[to_msat_solver.transfer_term(uf)] = uf;
  }

  // get mathsat terms for transition system
  msat_term msat_init = static_pointer_cast<MsatTerm>(
                            to_msat_solver.transfer_term(ts_.init(), BOOL))
                            ->get_msat_term();
  msat_term msat_trans = static_pointer_cast<MsatTerm>(
                             to_msat_solver.transfer_term(ts_.trans(), BOOL))
                             ->get_msat_term();
  msat_term msat_prop =
      static_pointer_cast<MsatTerm>(
          to_msat_solver.transfer_term(solver_->make_term(Not, bad_), BOOL))
          ->get_msat_term();
  unordered_map<msat_term, msat_term> msat_statevars;
  for (const auto & sv : ts_.statevars()) {
    msat_statevars[static_pointer_cast<MsatTerm>(
                       to_msat_solver.transfer_term(sv))
                       ->get_msat_term()] =
        static_pointer_cast<MsatTerm>(
            to_msat_solver.transfer_term(ts_.next(sv)))
            ->get_msat_term();
  }
  // initialize the transition system
  // NOTE: assuming not a liveprop
  ic3ia_ts.initialize(msat_statevars, msat_init, msat_trans, msat_prop, false);

  // just using default options for now
  ic3ia::Options ic3ia_opts;
  // the only options we pass through are
  // verbosity and random seed
  ic3ia_opts.seed = options_.random_seed_;
  ic3ia_opts.verbosity = options_.verbosity_;
  ic3ia::Logger & l = ic3ia::Logger::get();
  l.set_verbosity(ic3ia_opts.verbosity);
  // NOTE: assuming no LTL / liveness -- just adding because required
  ic3ia::LiveEncoder liveenc(ic3ia_ts, ic3ia_opts);
  ic3ia::IC3 ic3(ic3ia_ts, ic3ia_opts, liveenc);
  logger.log(1, "Running open-source ic3ia as backend.");
  msat_truth_value res = ic3.prove();

  if (res == MSAT_UNDEF) {
    return ProverResult::UNKNOWN;
  } else if (res == MSAT_TRUE) {
    invar_ = solver_->make_term(true);

    vector<ic3ia::TermList> ic3ia_invar;
    ic3.witness(ic3ia_invar);

    for (const auto & msat_clause : ic3ia_invar) {
      Term clause = solver_->make_term(false);
      assert(msat_clause.size());
      for (const msat_term & l : msat_clause) {
        clause = solver_->make_term(
            Or,
            clause,
            to_ts_solver.transfer_term(make_shared<MsatTerm>(env, l)));
      }
      invar_ = solver_->make_term(And, invar_, clause);
    }

    return ProverResult::TRUE;
  } else {
    assert(res == MSAT_FALSE);
    compute_witness(env, ic3, to_ts_solver);
    return ProverResult::FALSE;
  }
}

ProverResult MsatIC3IA::check_until(int k)
{
  throw PonoException("MsatIC3IA only supports prove, not check_until");
}

bool MsatIC3IA::compute_witness(msat_env env,
                                ic3ia::IC3 & ic3,
                                TermTranslator & to_ts_solver)
{
  // compute the witness, guided by the one from ic3ia
  // ic3ia does not give assignments to inputs, which we require
  assert(solver_->get_solver_enum() == MSAT);
  solver_->reset_assertions();

  Sort boolsort = solver_->make_sort(BOOL);

  vector<ic3ia::TermList> ic3ia_wit;
  ic3.witness(ic3ia_wit);
  assert(ic3ia_wit.size());

  // set reached_k_ so that it matches the counterexample length
  // reached_k_ was last bound without a counterexample, so
  // it's size - 2
  reached_k_ = ic3ia_wit.size() - 2;

  // set up a BMC query
  // with state variables constrained at each step
  solver_->assert_formula(unroller_.at_time(ts_.init(), 0));
  // assert that bad_ is not null
  // i.e. this prover was correctly initialized
  assert(bad_);
  solver_->assert_formula(unroller_.at_time(bad_, ic3ia_wit.size() - 1));
  for (size_t i = 0; i < ic3ia_wit.size(); ++i) {
    if (i + 1 < ic3ia_wit.size()) {
      solver_->assert_formula(unroller_.at_time(ts_.trans(), i));
    }

    for (const auto & msat_eq : ic3ia_wit[i]) {
      // create an smt-switch term for the equality
      Term eq = make_shared<MsatTerm>(env, msat_eq);
      Term solver_eq = to_ts_solver.transfer_term(eq);
      // either a literal or an equality
      assert(is_lit(solver_eq, boolsort) || solver_eq->get_op() == Equal);
      solver_->assert_formula(unroller_.at_time(solver_eq, i));
    }
  }

  Result r = solver_->check_sat();
  assert(r.is_sat());  // expecting a counterexample trace

  // rely on default compute_witness method to get model from solver_
  super::compute_witness();

  assert(ic3ia_wit.size() == witness_.size());
  return true;
}

}  // namespace pono
