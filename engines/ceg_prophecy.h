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
#include "refiners/array_axiom_enumerator.h"

namespace pono {

class CegProphecy : public CEGAR
{
  typedef CEGAR super;

  CegProphecy(const Property & p, smt::SolverEnum se);
  CegProphecy(const Property & p, const smt::SmtSolver & solver);
  CegProphecy(const PonoOptions & opt, const Property & p, smt::SolverEnum se);
  CegProphecy(const PonoOptions & opt,
              const Property & p,
              const smt::SmtSolver & solver);

 public:
  ProverResult check_until(int k) override;

 protected:
  const TransitionSystem & conc_ts_;
  const smt::SmtSolver & solver_;
  TransitionSystem abs_ts_;

  ArrayAbstractor aa_;
  ArrayAxiomEnumerator aae_;

  void abstract() override;
  void refine() override;
};

}  // namespace pono
