/*********************                                                  */
/*! \file ic3ia.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief IC3 via Implicit Predicate Abstraction (IC3IA) implementation. Based
*on
**        IC3 Modulo Theories via Implicit Predicate Abstraction
**            -- Alessandro Cimatti, Alberto Griggio, Sergio Mover, Stefano
*Tonetta
**        and the open source implementation:
**        https://es-static.fbk.eu/people/griggio/ic3ia/index.html
**/

#include "engines/mbic3.h"
#include "modifiers/implicit_predicate_abstractor.h"

namespace pono {

// process is mostly the same as model based IC3,
// but overrides a few methods to only use cubes/clauses over predicates
class IC3IA : public ModelBasedIC3
{
 public:
  IC3IA(Property & p, smt::SolverEnum se);
  IC3IA(Property & p, const smt::SmtSolver & slv);
  IC3IA(const PonoOptions & opt, Property & p, smt::SolverEnum se);
  IC3IA(const PonoOptions & opt, Property & p, const smt::SmtSolver & slv);
  virtual ~IC3IA();

  typedef ModelBasedIC3 super;

 protected:
  // TODO figure out which methods need to be overloaded here

  // TODO figure out relevant state

  RelationalTransitionSystem abs_ts_;

  ImplicitPredicateAbstractor ia_;

  smt::Term predabs_label_;
  ///< label that activates predicate equalities over next and abstract
  ///< statevars the same as ia_.predabs_label()

  smt::TermVec pred_statevars_;
  ///< state variables that IC3IA sees
  ///< these represent predicates and should stay aligned with the
  ///< predicate vector in the abstraction (ia_.predicates())
  ///< e.g. pred_statevars_[i] corresponds to ia_.predicates()[i]
};

}  // namespace pono
