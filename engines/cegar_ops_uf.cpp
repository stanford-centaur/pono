/*********************                                                        */
/*! \file cegar_ops_uf.cpp
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief A simple CEGAR loop that abstracts operators with uninterpreted
**        functions and refines by concretizing the operators
**
**/

#include "engines/cegar_ops_uf.h"

#include "core/rts.h"
#include "core/ts.h"
#include "engines/bmc.h"
#include "engines/bmc_simplepath.h"
#include "engines/ceg_prophecy_arrays.h"
#include "engines/dual_approx_reach.h"
#include "engines/ic3ia.h"
#include "engines/ic3sa.h"
#include "engines/interp_seq_mc.h"
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"
#include "smt/available_solvers.h"
#include "utils/exceptions.h"
#include "utils/logger.h"
#include "utils/ts_manipulation.h"

using namespace smt;
using namespace std;

namespace pono {

template <class Prover_T>
CegarOpsUf<Prover_T>::CegarOpsUf(const SafetyProperty & p,
                                 const TransitionSystem & ts,
                                 const SmtSolver & solver,
                                 PonoOptions opt)
    : super(p, create_fresh_ts(ts.is_functional(), solver), solver, opt),
      conc_ts_(ts, super::to_prover_solver_),
      prover_ts_(super::prover_interface_ts()),
      oa_(conc_ts_, prover_ts_, opt.ceg_bv_arith_as_free_symbol_),
      cegopsuf_solver_(
          create_solver(solver->get_solver_enum(), opt.logging_smt_solver_)),
      to_cegopsuf_solver_(cegopsuf_solver_),
      from_cegopsuf_solver_(super::prover_interface_ts().solver()),
      cegopsuf_ts_(cegopsuf_solver_),
      cegopsuf_un_(cegopsuf_ts_)
{
  cegopsuf_solver_->set_opt("produce-unsat-assumptions", "true");

  // store the original ts for later use in witness generation
  super::orig_ts_ = ts;
}

template <class Prover_T>
void CegarOpsUf<Prover_T>::set_ops_to_abstract(
    const UnorderedOpSet & ops_to_abstract)
{
  oa_.set_ops_to_abstract(ops_to_abstract);
}

template <class Prover_T>
ProverResult CegarOpsUf<Prover_T>::check_until(int k)
{
  initialize();

  ProverResult res = ProverResult::FALSE;
  while (res == ProverResult::FALSE) {
    // need to call parent's check_until in case it
    // is another cegar loop rather than an engine
    res = super::check_until(k);

    if (res == ProverResult::FALSE) {
      if (!cegar_refine()) {
        return ProverResult::FALSE;
      }
    }
  }

  if (res == ProverResult::TRUE && super::invar_) {
    if (super::options_.check_invar_) {
      // TODO: Find an example that triggers the error and
      // document it in the issue tracker.
      logger.log(1,
                 "WARNING: Concretizing invariant in CegarOpsUf... "
                 "This feature is not stable and may encounter errors.");
      super::invar_ = oa_.concrete(super::invar_);
    } else {
      super::invar_ = NULL;
    }
  }

  return res;
}

template <class Prover_T>
void CegarOpsUf<Prover_T>::initialize()
{
  if (super::initialized_) {
    return;
  }

  // specify which cegar_abstract in case
  // we're inheriting from another cegar algorithm
  CegarOpsUf::cegar_abstract();
  super::initialize();

  // update local version of ts over fresh solver
  cegopsuf_ts_ = TransitionSystem(prover_ts_, to_cegopsuf_solver_);
  // update cache
  UnorderedTermMap & cache = from_cegopsuf_solver_.get_cache();
  Term nv;
  for (const auto & sv : prover_ts_.statevars()) {
    nv = prover_ts_.next(sv);
    cache[to_cegopsuf_solver_.transfer_term(sv)] = sv;
    cache[to_cegopsuf_solver_.transfer_term(nv)] = nv;
  }

  for (const auto & iv : prover_ts_.inputvars()) {
    cache[to_cegopsuf_solver_.transfer_term(iv)] = iv;
  }

  const UnorderedTermMap & abs_terms = oa_.abstract_terms();
  Sort boolsort = cegopsuf_solver_->make_sort(BOOL);
  Term lbl;
  for (const auto & elem : abs_terms) {
    lbl = cegopsuf_solver_->make_symbol(
        "cegopsuf_assump_" + std::to_string(elem.second->get_id()), boolsort);
    cegopsuf_labels_[to_cegopsuf_solver_.transfer_term(elem.first)] = lbl;
    cache[to_cegopsuf_solver_.transfer_term(elem.first)] = elem.first;
    cache[to_cegopsuf_solver_.transfer_term(elem.second)] = elem.second;
  }
}

template <class Prover_T>
void CegarOpsUf<Prover_T>::cegar_abstract()
{
  oa_.do_abstraction();

  assert(super::bad_);
  super::bad_ = oa_.abstract(super::bad_);
  cegopsuf_bad_ = to_cegopsuf_solver_.transfer_term(super::bad_, BOOL);
}

template <class Prover_T>
bool CegarOpsUf<Prover_T>::cegar_refine()
{
  const UnorderedTermMap & abs_terms = oa_.abstract_terms();
  if (abs_terms.size() == 0) {
    // no abstraction, let prover compute witness
    Prover_T::compute_witness();
    return false;
  }

  size_t cex_length = super::witness_length();

  // create bmc formula for abstract system
  Term bmcform = cegopsuf_un_.at_time(cegopsuf_ts_.init(), 0);
  for (size_t i = 0; i < cex_length; ++i) {
    bmcform = cegopsuf_solver_->make_term(
        And, bmcform, cegopsuf_un_.at_time(cegopsuf_ts_.trans(), i));
  }
  bmcform = cegopsuf_solver_->make_term(
      And, bmcform, cegopsuf_un_.at_time(cegopsuf_bad_, cex_length));

  cegopsuf_solver_->push();
  cegopsuf_solver_->assert_formula(bmcform);

  TermVec assumps;
  TermVec equalities;
  for (const auto & elem : abs_terms) {
    Term l = to_cegopsuf_solver_.transfer_term(elem.first);
    Term r = to_cegopsuf_solver_.transfer_term(elem.second);

    assert(cegopsuf_labels_.find(l) != cegopsuf_labels_.end());
    Term lbl = cegopsuf_labels_.at(l);
    Term uf_eq = cegopsuf_solver_->make_term(Equal, l, r);

    Term unrolled_uf_eq = cegopsuf_un_.at_time(uf_eq, 0);
    for (size_t i = 1; i < cex_length; ++i) {
      unrolled_uf_eq = cegopsuf_solver_->make_term(
          And, unrolled_uf_eq, cegopsuf_un_.at_time(uf_eq, i));
    }
    if (cegopsuf_ts_.only_curr(uf_eq)) {
      unrolled_uf_eq = cegopsuf_solver_->make_term(
          And, unrolled_uf_eq, cegopsuf_un_.at_time(uf_eq, cex_length));
    }

    equalities.push_back(uf_eq);
    Term imp = cegopsuf_solver_->make_term(Implies, lbl, unrolled_uf_eq);
    cegopsuf_solver_->assert_formula(imp);
    assumps.push_back(lbl);
  }
  assert(assumps.size() == equalities.size());

  Result r = cegopsuf_solver_->check_sat_assuming(assumps);

  // do refinement if needed
  if (r.is_unsat()) {
    UnorderedTermSet axioms;
    UnorderedTermSet core;
    cegopsuf_solver_->get_unsat_assumptions(core);

    for (size_t i = 0; i < assumps.size(); ++i) {
      if (core.find(assumps[i]) != core.end()) {
        Term eq = equalities[i];
        axioms.insert(eq);
        logger.log(2, "CegarOpsUf adding refinement axiom {}", eq);
        // need to refine both systems
        if (cegopsuf_ts_.is_functional()) {
          cegopsuf_ts_.add_constraint(eq);
        } else {
          RelationalTransitionSystem & ts =
              static_cast<RelationalTransitionSystem &>(cegopsuf_ts_);
          if (ts.only_curr(eq) && cex_length == 0) {
            ts.constrain_init(eq);
          }

          ts.constrain_trans(eq);
          if (ts.no_next(eq)) {
            ts.constrain_trans(ts.next(eq));
          }
        }
      }
    }
    refine_subprover_ts(axioms, cex_length > 0);
  } else if (r.is_sat()) {
    // store CEX trace
    store_witness();
  } else {
    throw PonoException("CegarOpsUf::cegar_refine: solver returned "
                        + r.to_string());
  }
  cegopsuf_solver_->pop();

  return r.is_unsat();
}

template <class Prover_T>
void CegarOpsUf<Prover_T>::store_witness()
{
  auto & witness = super::witness_;
  witness.clear();
  witness.reserve(super::witness_length() + 1);

  for (size_t i = 0, wl = super::witness_length(); i <= wl; ++i) {
    witness.push_back(UnorderedTermMap());
    UnorderedTermMap & map = witness.back();

    for (const auto & v : cegopsuf_ts_.statevars()) {
      const Term & vi = cegopsuf_un_.at_time(v, i);
      const Term & r = cegopsuf_solver_->get_value(vi);
      map[from_cegopsuf_solver_.transfer_term(v)] =
          from_cegopsuf_solver_.transfer_term(r);
    }

    for (const auto & v : cegopsuf_ts_.inputvars()) {
      const Term & vi = cegopsuf_un_.at_time(v, i);
      const Term & r = cegopsuf_solver_->get_value(vi);
      map[from_cegopsuf_solver_.transfer_term(v)] =
          from_cegopsuf_solver_.transfer_term(r);
    }

    for (const auto & elem : cegopsuf_ts_.named_terms()) {
      const Term & ti = cegopsuf_un_.at_time(elem.second, i);
      const Term & r = cegopsuf_solver_->get_value(ti);
      map[from_cegopsuf_solver_.transfer_term(elem.second)] =
          from_cegopsuf_solver_.transfer_term(r);
    }
  }
}

template <class Prover_T>
void CegarOpsUf<Prover_T>::refine_subprover_ts_base(
    const UnorderedTermSet & axioms, bool skip_init)
{
  for (const auto & a : axioms) {
    Term ta = from_cegopsuf_solver_.transfer_term(a, BOOL);

    if (prover_ts_.is_functional()) {
      prover_ts_.add_constraint(ta);
    } else {
      RelationalTransitionSystem & ts =
          static_cast<RelationalTransitionSystem &>(prover_ts_);
      if (ts.only_curr(ta) && !skip_init) {
        ts.constrain_init(ta);
      }

      ts.constrain_trans(ta);
      if (ts.no_next(ta)) {
        ts.constrain_trans(ts.next(ta));
      }
    }
  }
}

template <class Prover_T>
void CegarOpsUf<Prover_T>::refine_subprover_ts(const UnorderedTermSet & axioms,
                                               bool skip_init)
{
  refine_subprover_ts_base(axioms, skip_init);
  super::reset_env();
}

template <>
void CegarOpsUf<IC3IA>::refine_subprover_ts(const UnorderedTermSet & axioms,
                                            bool skip_init)
{
  refine_subprover_ts_base(axioms, skip_init);
  super::reabstract();
}

template <>
void CegarOpsUf<IC3SA>::refine_subprover_ts(const UnorderedTermSet & axioms,
                                            bool skip_init)
{
  for (const auto & a : axioms) {
    Term ta = from_cegopsuf_solver_.transfer_term(a, BOOL);

    // mine for new terms for the term abstraction
    super::add_to_term_abstraction(ta);

    assert(prover_ts_.is_functional());
    prover_ts_.add_constraint(ta);

    // main ts_ in IC3SA is a relational view of the original
    // (functional) ts
    // for now, easiest to just add it directly to that system
    // as well
    // TODO clean this up later
    RelationalTransitionSystem & ts =
        static_cast<RelationalTransitionSystem &>(super::ts_);
    if (ts.only_curr(ta) && !skip_init) {
      ts.constrain_init(ta);
      super::solver_->assert_formula(
          solver_->make_term(Implies, super::init_label_, ta));
    }

    ts.constrain_trans(ta);
    super::solver_->assert_formula(
        solver_->make_term(Implies, super::trans_label_, ta));

    if (ts.no_next(ta)) {
      // NOTE: don't need to add next version to term abstraction
      // because only keeps current state vars anyway
      Term next_ta = ts.next(ta);
      ts.constrain_trans(next_ta);
      super::solver_->assert_formula(
          solver_->make_term(Implies, super::trans_label_, next_ta));
    }
  }
}

template <>
void CegarOpsUf<CegProphecyArrays<IC3IA>>::refine_subprover_ts(
    const UnorderedTermSet & axioms, bool skip_init)
{
  super::refine_ts(axioms);
}

// TODO add other template classes
template class CegarOpsUf<Bmc>;
template class CegarOpsUf<BmcSimplePath>;
template class CegarOpsUf<DualApproxReach>;
template class CegarOpsUf<IC3IA>;
template class CegarOpsUf<IC3SA>;
template class CegarOpsUf<InterpolantMC>;
template class CegarOpsUf<InterpSeqMC>;
template class CegarOpsUf<KInduction>;
template class CegarOpsUf<CegProphecyArrays<IC3IA>>;

}  // namespace pono
