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

using namespace smt;
using namespace std;

namespace pono {

template <class Prover_T>
CegarOpsUf<Prover_T>::CegarOpsUf(const Property & p,
                                 const TransitionSystem & ts,
                                 const SmtSolver & solver,
                                 PonoOptions opt)
  : super(p, create_fresh_ts(ts.is_functional(), solver), solver, opt),
    conc_ts_(ts, super::to_prover_solver_),
    prover_ts_(super::prover_interface_ts()),
    oa_(conc_ts_, prover_ts_),
    cegopsuf_solver_(create_solver(solver->get_solver_enum())),
    to_cegopsuf_solver_(cegopsuf_solver_),
    from_cegopsuf_solver_(super::prover_interface_ts().solver()),
    cegopsuf_ts_(cegopsuf_solver_),
    cegopsuf_un_(cegopsuf_ts_)
{
  cegopsuf_solver_->set_opt("prove-unsat-cores", "true");
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
    lbl = cegopsuf_solver_->make_symbol("__cegopsuf_assump_" + elem.second->to_string(),
                                        boolsort);
    cegopsuf_labels_[elem.first] = lbl;
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
    Term lbl = cegopsuf_labels_.at(elem.first);
    assumps.push_back(lbl);
    Term uf_eq = cegopsuf_solver_->make_term(Equal, elem.first, elem.second);

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
  }
  assert(assumps.size() == equalities.size());

  Result r = cegopsuf_solver_->check_sat_assuming(assumps);

  // do refinement if needed
  if (r.is_unsat()) {
    UnorderedTermSet core;
    cegopsuf_solver_->get_unsat_core(core);
    for (size_t i = 0; i < assumps.size(); ++i) {
      if (core.find(assumps[i]) != core.end()) {
        logger.log(2, "CegarValues adding refinement axiom {}", equalities[i]);
        // need to refine both systems
        cegopsuf_ts_.add_constraint(equalities[i]);
        // TODO this should be more modular
        //      can't assume super::abs_ts_ is the right one to constrain
        refine_subprover_ts(equalities[i], cex_length > 0);
      }
    }
  }
  cegopsuf_solver_->pop();

  return r.is_unsat();
}

template <class Prover_T>
void CegarOpsUf<Prover_T>::refine_subprover_ts(const Term & constraint,
                                               bool skip_init)
{
  throw PonoException("CegarValues::refine_subprover_ts NYI for generic case");
}

template <>
void CegarOpsUf<IC3IA>::refine_subprover_ts(const Term & constraint,
                                            bool skip_init)
{
  Term transferred_constraint =
    from_cegopsuf_solver_.transfer_term(constraint, BOOL);

  assert(!ts_.is_functional());
  RelationalTransitionSystem & rts =
    static_cast<RelationalTransitionSystem &>(ts_);

  rts.constrain_trans(transferred_constraint);
  if (rts.no_next(transferred_constraint)) {
    rts.constrain_trans(rts.next(transferred_constraint));
  }

  if (rts.only_curr(transferred_constraint) && !skip_init) {
    rts.constrain_init(transferred_constraint);
  }
}

// TODO add other template classes
template class CegarOpsUf<IC3IA>;

} // namespace pono
