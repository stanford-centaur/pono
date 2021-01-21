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

#include "core/fts.h"
#include "core/rts.h"

#include "engines/ic3ia.h"

#include "smt/available_solvers.h"

#include "utils/exceptions.h"
#include "utils/logger.h"
#include "utils/make_provers.h"
#include "utils/ts_manipulation.h"
#include "utils/term_analysis.h"

using namespace smt;
using namespace std;

namespace pono {

template <class Prover_T>
CegarOpsUf<Prover_T>::CegarOpsUf(const Property & p,
                                 const TransitionSystem & ts,
                                 const SmtSolver & solver,
                                 PonoOptions opt)
  : super(p, create_fresh_ts(false, solver), solver, opt),
    conc_ts_(ts, super::to_prover_solver_),
    prover_ts_(super::prover_interface_ts()),
    oa_(conc_ts_, prover_ts_),
    cegopsuf_solver_(create_solver(solver->get_solver_enum())),
    to_cegopsuf_solver_(cegopsuf_solver_),
    from_cegopsuf_solver_(super::prover_interface_ts().solver()),
    cegopsuf_ts_(cegopsuf_solver_),
    cegopsuf_un_(cegopsuf_ts_)
{
  cegopsuf_solver_->set_opt("produce-unsat-cores", "true");
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
    // update the invariant

    //TODO
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
    lbl = cegopsuf_solver_->make_symbol("cegopsuf_assump_" + std::to_string(elem.second->get_id()),
                                        boolsort);
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

  const UnorderedTermMap & abs_terms = oa_.abstract_terms();
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
      unrolled_uf_eq = cegopsuf_solver_->make_term(And, unrolled_uf_eq,
                                                   cegopsuf_un_.at_time(uf_eq, i));
    }
    if (cegopsuf_ts_.only_curr(uf_eq)) {
      unrolled_uf_eq = cegopsuf_solver_->make_term(And, unrolled_uf_eq,
                                                   cegopsuf_un_.at_time(uf_eq, cex_length));
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
    cegopsuf_solver_->get_unsat_core(core);

    assert(!cegopsuf_ts_.is_functional());
    RelationalTransitionSystem & ts =
      static_cast<RelationalTransitionSystem &>(cegopsuf_ts_);
    for (size_t i = 0; i < assumps.size(); ++i) {
      if (core.find(assumps[i]) != core.end()) {
        Term eq = equalities[i];
        axioms.insert(eq);
        logger.log(2, "CegarOpsUf adding refinement axiom {}", eq);
        // need to refine both systems
        ts.constrain_trans(eq);
        if (ts.no_next(eq)) {
          ts.constrain_trans(ts.next(eq));
        }
        if (cex_length == 0 && ts.only_curr(eq)) {
          ts.constrain_init(eq);
        }
      }
    }
    refine_subprover_ts(axioms, cex_length > 0);
  }
  cegopsuf_solver_->pop();

  return r.is_unsat();
}

template <class Prover_T>
void CegarOpsUf<Prover_T>::refine_subprover_ts(const UnorderedTermSet & axioms,
                                               bool skip_init)
{
  throw PonoException("CegarOpsUf::refine_subprover_ts NYI for generic case");
}

template <>
void CegarOpsUf<IC3IA>::refine_subprover_ts(const UnorderedTermSet & axioms,
                                            bool skip_init)
{
  assert(!prover_ts_.is_functional());
  RelationalTransitionSystem & rts =
    static_cast<RelationalTransitionSystem &>(prover_ts_);

  for (const auto & a : axioms) {
    Term ta = from_cegopsuf_solver_.transfer_term(a, BOOL);

    rts.constrain_trans(ta);
    if (rts.no_next(ta)) {
      rts.constrain_trans(rts.next(ta));
    }

    if (rts.only_curr(ta) && !skip_init) {
      rts.constrain_init(ta);
    }
  }

  // add predicates from init and bad
  UnorderedTermSet preds;
  get_predicates(super::solver_, super::ts_.init(), preds, false, false, true);
  get_predicates(super::solver_, super::bad_, preds, false, false, true);
  get_predicates(super::solver_, super::get_frame_term(1), preds, false, false, true);
  super::predset_.clear();
  super::predlbls_.clear();

  // don't add boolean symbols that are never used in the system
  // this is an optimization and a fix for some options
  // if using mathsat with bool_model_generation
  // it will fail to get the value of symbols that don't
  // appear in the query
  // thus we don't include those symbols in our cubes
  UnorderedTermSet used_symbols;
  get_free_symbolic_consts(super::ts_.init(), used_symbols);
  get_free_symbolic_consts(super::ts_.trans(), used_symbols);
  get_free_symbolic_consts(super::bad_, used_symbols);

  // reset init and trans -- done with calling ia_.do_abstraction
  // then add all boolean constants as (precise) predicates
  for (const auto & p : super::ia_.do_abstraction()) {
    preds.insert(p);
  }

  // reset the solver
  super::reset_solver();

  // add predicates
  for (const auto &p : preds) {
    super::add_predicate(p);
  }
}

// TODO add other template classes
template class CegarOpsUf<IC3IA>;

} // namespace pono
