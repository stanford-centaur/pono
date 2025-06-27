/*********************                                                        */
/*! \file cegar.h
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief An abstract base class for Counter-Example Guided Abstraction
*Refinement
**        (CEGAR) techniques.
**
**/

#pragma once

#include "core/prop.h"
#include "core/ts.h"
#include "options/options.h"
#include "smt-switch/smt.h"

namespace pono {

template <class Prover_T>
class CEGAR : public Prover_T
{
  typedef Prover_T super;

 public:
  CEGAR(const SafetyProperty & p,
        const TransitionSystem & ts,
        const smt::SmtSolver & solver,
        PonoOptions opt = PonoOptions())
      : super(p, ts, solver, opt)
  {
  }

  virtual ~CEGAR() {}

 protected:
  /** Abstract the transition system -- usually only performed once
   *  (in initialize)
   */
  virtual void cegar_abstract() = 0;
  /** Refine the abstracted transition system
   *  Typically performed in a refinement loop
   *  @return true iff it was successfully refined
   */
  virtual bool cegar_refine() = 0;
};

}  // namespace pono
