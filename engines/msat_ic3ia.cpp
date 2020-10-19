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

#include "assert.h"
#include "smt-switch/msat_solver.h"
#include "utils/logger.h"

using namespace ic3ia;
using namespace smt;
using namespace std;

namespace pono {

MsatIC3IA::MsatIC3IA(Property & p, smt::SolverEnum se) : super(p, se) {}

MsatIC3IA::MsatIC3IA(Property & p, const SmtSolver & solver) : super(p, solver)
{
}

MsatIC3IA::MsatIC3IA(const PonoOptions & opt, Property & p, smt::SolverEnum se)
  : super(opt, p, se) {}

MsatIC3IA::MsatIC3IA(const PonoOptions & opt, Property & p, const SmtSolver & solver)
  : super(opt, p, solver)
{
}

ProverResult MsatIC3IA::prove()
{
  if (ts_.solver()->get_solver_enum() != MSAT) {
    throw PonoException("MsatIC3IA only supports mathsat solver.");
  }
  const SmtSolver & solver = ts_.solver();
  shared_ptr<const MsatSolver> msat_solver =
      static_pointer_cast<MsatSolver>(solver);
  ::ic3ia::TransitionSystem ic3ia_ts(msat_solver->get_msat_env());

  // get mathsat terms for transition system
  msat_term msat_init =
      static_pointer_cast<MsatTerm>(ts_.init())->get_msat_term();
  msat_term msat_trans =
      static_pointer_cast<MsatTerm>(ts_.trans())->get_msat_term();
  msat_term msat_prop =
      static_pointer_cast<MsatTerm>(property_.prop())->get_msat_term();
  unordered_map<msat_term, msat_term> msat_statevars;
  for (auto sv : ts_.statevars()) {
    msat_statevars[static_pointer_cast<MsatTerm>(sv)->get_msat_term()] =
        static_pointer_cast<MsatTerm>(ts_.next(sv))->get_msat_term();
  }
  // initialize the transition system
  // NOTE: assuming not a liveprop
  ic3ia_ts.initialize(msat_statevars, msat_init, msat_trans, msat_prop, false);

  // just using default options for now
  ic3ia::Options ic3ia_opts;
  // the only option we pass through is verbosity
  ic3ia_opts.verbosity = options_.verbosity_;
  ic3ia::Logger & l = ic3ia::Logger::get();
  l.set_verbosity(ic3ia_opts.verbosity);
  // NOTE: assuming no LTL / liveness -- just adding because required
  LiveEncoder liveenc(ic3ia_ts, ic3ia_opts);
  IC3 ic3(ic3ia_ts, ic3ia_opts, liveenc);
  logger.log(1, "Running open-source ic3ia as backend.");
  msat_truth_value res = ic3.prove();

  if (res == MSAT_UNDEF) {
    return ProverResult::UNKNOWN;
  } else if (res == MSAT_TRUE) {
    // TODO populate invar_
    return ProverResult::TRUE;
  } else {
    assert(res == MSAT_FALSE);
    compute_witness(ic3);
    return ProverResult::FALSE;
  }
}

ProverResult MsatIC3IA::check_until(int k)
{
  throw PonoException("MsatIC3IA only supports prove, not check_until");
}

bool MsatIC3IA::compute_witness(ic3ia::IC3 & ic3)
{
  vector<ic3ia::TermList> ic3ia_wit;
  ic3.witness(ic3ia_wit);

  const SmtSolver & solver = ts_.solver();
  shared_ptr<const MsatSolver> msat_solver =
      static_pointer_cast<MsatSolver>(solver);
  msat_env env = msat_solver->get_msat_env();

  for (auto terms : ic3ia_wit) {
    witness_.push_back(UnorderedTermMap());
    UnorderedTermMap & map = witness_.back();

    for (auto msat_eq : terms) {
      // create an smt-switch term
      Term eq = make_shared<MsatTerm>(env, msat_eq);
      TermVec children(eq->begin(), eq->end());
      assert(children.size() == 2);
      Term sym = children[0];
      Term val = children[1];

      if (!sym->is_symbolic_const()) {
        // got the order wrong, reverse it
        std::swap(sym, val);
      }
      map[sym] = val;
    }
  }

  assert(ic3ia_wit.size() == witness_.size());
  return true;
}

}  // namespace pono
