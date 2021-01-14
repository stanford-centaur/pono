/*********************                                                        */
/*! \file cegar_values.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief A simple CEGAR loop that abstracts values with frozen variables
**        and refines by constraining the variable to the value again
**
**/

#include "engines/cegar_values.h"

#include <unordered_set>

#include "core/fts.h"
#include "core/rts.h"
#include "engines/ceg_prophecy_arrays.h"
#include "engines/ic3ia.h"
#include "smt-switch/identity_walker.h"
#include "smt/available_solvers.h"
#include "utils/exceptions.h"
#include "utils/make_provers.h"

using namespace smt;
using namespace std;

namespace pono {

// TODO add a value abstractor
//      make sure not to introduce nonlinearities
//      implement generic backend
//      then specialize for IC3IA

static unordered_set<PrimOp> nl_ops({ Mult,
                                      Div,
                                      Mod,
                                      Abs,
                                      Pow,
                                      IntDiv,
                                      BVMul,
                                      BVUdiv,
                                      BVSdiv,
                                      BVUrem,
                                      BVSrem,
                                      BVSmod });

/** Class for abstracting values with frozen variables
 */
class ValueAbstractor : public smt::IdentityWalker
{
 public:
  ValueAbstractor(TransitionSystem & ts, UnorderedTermMap & abstracted_values)
      : smt::IdentityWalker(ts.solver(), false),
        ts_(ts),
        abstracted_values_(abstracted_values)
  {
  }

 protected:
  smt::WalkerStepResult visit_term(smt::Term & term) override
  {
    if (!preorder_) {
      if (term->is_value() && term->get_sort()->get_sort_kind() != ARRAY) {
        // create a frozen variable
        Term frozen_var =
            ts_.make_statevar("__abs_" + term->to_string(), term->get_sort());
        // make sure it's frozen
        ts_.assign_next(frozen_var, frozen_var);
        save_in_cache(term, frozen_var);
        abstracted_values_[frozen_var] = term;
        return Walker_Continue;
      }

      Op op = term->get_op();
      if (!op.is_null() && nl_ops.find(op.prim_op) == nl_ops.end()) {
        // only rebuild terms with operators that can't
        // create nonlinearities by replacing constant values
        TermVec cached_children;
        Term c;
        for (const auto & t : term) {
          c = t;
          query_cache(t, c);
          cached_children.push_back(c);
        }
        save_in_cache(term, solver_->make_term(op, cached_children));
      } else {
        save_in_cache(term, term);
      }
    }
    return Walker_Continue;
  }

  TransitionSystem & ts_;
  UnorderedTermMap & abstracted_values_;
};

TransitionSystem create_fresh_ts(bool functional, const SmtSolver & solver)
{
  if (functional) {
    return FunctionalTransitionSystem(solver);
  } else {
    return RelationalTransitionSystem(solver);
  }
}

template <class Prover_T>
CegarValues<Prover_T>::CegarValues(const Property & p,
                                   const TransitionSystem & ts,
                                   const smt::SmtSolver & solver,
                                   PonoOptions opt)
    : super(p, create_fresh_ts(ts.is_functional(), solver), solver, opt),
      conc_ts_(ts, super::to_prover_solver_),
      prover_ts_(super::prover_interface_ts())
{
}

template <class Prover_T>
ProverResult CegarValues<Prover_T>::check_until(int k)
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

  return res;
}

template <class Prover_T>
void CegarValues<Prover_T>::initialize()
{
  if (super::initialized_) {
    return;
  }

  // specify which cegar_abstract in case
  // we're inheriting from another cegar algorithm
  CegarValues::cegar_abstract();
  super::initialize();
}

template <class Prover_T>
void CegarValues<Prover_T>::cegar_abstract()
{
  prover_ts_ = conc_ts_;
  ValueAbstractor va(prover_ts_, to_vals_);

  // now update with abstraction
  if (prover_ts_.is_functional()) {
    throw PonoException("Functional TS NYI in cegar_values");
  } else {
    Term init = prover_ts_.init();
    Term trans = prover_ts_.trans();
    static_cast<RelationalTransitionSystem &>(prover_ts_)
        .set_behavior(va.visit(init), va.visit(trans));
  }

  // TODO clean this up
  // kind of complicated to get property because super initialize
  // hasn't been run yet
  Term prop_term = (prover_ts_.solver() == super::orig_property_.solver())
                       ? super::orig_property_.prop()
                       : super::to_prover_solver_.transfer_term(
                           super::orig_property_.prop(), BOOL);
  super::bad_ =
      super::solver_->make_term(smt::PrimOp::Not, va.visit(prop_term));

  assert(super::bad_);

  // expecting to have had values to abstract
  assert(to_vals_.size());
}

template <class Prover_T>
bool CegarValues<Prover_T>::cegar_refine()
{
  throw PonoException("NYI");
}

// TODO add other template classes
template class CegarValues<CegProphecyArrays<IC3IA>>;

}  // namespace pono
