/*********************                                                        */
/*! \file fault_injector.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the cosa2 project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Modifies a transition system such that faults can occur on state
**        variables at each timestep. This assumes a functional model, and
**        inserts a mux and a boolean FAULT variable for each state element.
**        If the fault variable is true, it can get a non-deterministic value
**
**/

#pragma once

#include "smt-switch/smt.h"

#include "fts.h"

namespace cosa
{

  class FaultInjector
  {
  public:
    FaultInjector(FunctionalTransitionSystem& fts) : fts_(fts), faulty_fts_(fts.solver())
    {
      do_fault_injection();
    };

    FunctionalTransitionSystem faulty_transition_system() const { return faulty_fts_; };

    smt::TermVec fault_sigs() const { return fault_sigs_; };

  protected:

    void do_fault_injection();

    FunctionalTransitionSystem & fts_;
    FunctionalTransitionSystem faulty_fts_;

    smt::UnorderedTermMap state2faultsel_;
    smt::UnorderedTermMap faultsel2state_;

    smt::UnorderedTermMap state2faultval_;
    smt::UnorderedTermMap faultval2state_;

    smt::UnorderedTermMap state2faultsig_;
    smt::UnorderedTermMap faultsig2state_;

    smt::TermVec fault_sigs_;
  };

}
