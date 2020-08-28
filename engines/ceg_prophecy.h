/*********************                                                        */
/*! \file ceg_prophecy.h
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

#pragma once

#include "engines/cegar.h"
#include "modifiers/array_abstractor.h"
#include "modifiers/prophecy_modifier.h"
#include "options/options.h"
#include "refiners/array_axiom_enumerator.h"

namespace pono {

class CegProphecy : public CEGAR
{
  typedef CEGAR super;

 public:
  CegProphecy(const Property & p, Engine e, smt::SolverEnum se);
  CegProphecy(const Property & p, Engine e, const smt::SmtSolver & solver);
  CegProphecy(const PonoOptions & opt,
              const Property & p,
              Engine e,
              smt::SolverEnum se);
  CegProphecy(const PonoOptions & opt,
              const Property & p,
              Engine e,
              const smt::SmtSolver & solver);

  ProverResult prove() override;
  ProverResult check_until(int k) override;

  void initialize() override;

 protected:
  const TransitionSystem & conc_ts_;
  const smt::SmtSolver & solver_;
  RelationalTransitionSystem abs_ts_;
  Engine e_;

  // TODO: see if there's a better organization where we can re-use the same
  // unroller currently this is very important, or it won't unroll the
  // abstracted variables correctly and this will fail SILENTLY
  Unroller abs_unroller_;  ///< unroller for abs_ts_
  ArrayAbstractor aa_;
  ArrayAxiomEnumerator aae_;
  ProphecyModifier pm_;

  size_t num_added_axioms_ =
      0;  ///< set by refine to the number of added axioms

  void abstract() override;
  bool refine() override;
};

}  // namespace pono
