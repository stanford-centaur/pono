/*********************                                                        */
/*! \file ceg_prophecy.cpp
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief An implementation of Counter-Example Guided Prophecy for array
**        model checking. It is parameterized by an underlying model checking
**        procedure which need not handle arrays (only UF). However, a common
**        instantiation is with an IC3-style procedure, in which case we
**        often refer to this algorithm as "prophic3".
**
**/

#include "assert.h"

#include "engines/ceg_prophecy.h"
#include "utils/make_provers.h"

using namespace smt;
using namespace std;

namespace pono {

CegProphecy::CegProphecy(const Property & p, Engine e, smt::SolverEnum se)
    : super(p, se),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      e_(e),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, unroller_)
{
  super::initialize();
}

CegProphecy::CegProphecy(const Property & p, Engine e, const SmtSolver & solver)
    : super(p, solver),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      e_(e),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, unroller_)
{
  super::initialize();
}

CegProphecy::CegProphecy(const PonoOptions & opt,
                         const Property & p,
                         Engine e,
                         smt::SolverEnum se)
    : super(opt, p, se),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      e_(e),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, unroller_)
{
  super::initialize();
}

CegProphecy::CegProphecy(const PonoOptions & opt,
                         const Property & p,
                         Engine e,
                         const smt::SmtSolver & solver)
    : super(opt, p, solver),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      e_(e),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, unroller_)
{
  super::initialize();
}

ProverResult CegProphecy::prove()
{
  ProverResult res = ProverResult::FALSE;
  while (res == ProverResult::FALSE) {
    // Refine the system
    // heuristic -- stop refining when no new axioms are needed.
    do {
      if (!refine()) {
        return ProverResult::FALSE;
      }
      reached_k_++;
    } while (num_added_axioms_);

    Property latest_prop(abs_ts_, solver_->make_term(Not, bad_));
    shared_ptr<Prover> prover =
        make_prover(e_, latest_prop, solver_->get_solver_enum(), options_);
    res = prover->prove();
  }

  return res;
}

ProverResult CegProphecy::check_until(int k)
{
  ProverResult res = ProverResult::FALSE;
  while (res == ProverResult::FALSE && reached_k_ <= k) {
    // Refine the system
    // heuristic -- stop refining when no new axioms are needed.
    do {
      if (!refine()) {
        return ProverResult::FALSE;
      }
      reached_k_++;
    } while (num_added_axioms_ && reached_k_ <= k);

    Property latest_prop(abs_ts_, solver_->make_term(Not, bad_));
    shared_ptr<Prover> prover =
        make_prover(e_, latest_prop, solver_->get_solver_enum(), options_);
    res = prover->check_until(k);
  }
  return res;
}

void CegProphecy::abstract()
{
  // this is a No-Op
  // the ArrayAbstractor already abstracted the transition system on
  // construction

  // the abstract system should have all the same state and input variables
  // but abstracted
  // plus it will have some new variables for the witnesses and lambdas
  assert(abs_ts_.statevars().size() >= conc_ts_.statevars().size());
  assert(abs_ts_.inputvars().size() >= conc_ts_.inputvars().size());
}

bool CegProphecy::refine()
{
  num_added_axioms_ = 0;
  // TODO use ArrayAxiomEnumerator and modifiers to refine the system

  throw PonoException("NYI");
}

}  // namespace pono
