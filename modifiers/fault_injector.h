/*********************                                                        */
/*! \file fault_injector.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
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

#include "utils/exceptions.h"
#include "core/fts.h"

namespace pono {

class FaultInjector
{
 public:
 FaultInjector(TransitionSystem & fts, smt::UnorderedTermSet to_ignore={})
   : fts_(fts), faulty_fts_(fts.solver()), statevars_to_ignore_(to_ignore)
  {
    if (!fts.is_functional()) {
      throw PonoException("Can only do fault injection on functional systems.");
    }
    do_fault_injection();
  };

  virtual ~FaultInjector(){};

  FunctionalTransitionSystem & faulty_transition_system()
  {
    return faulty_fts_;
  };

  smt::TermVec fault_sigs() const { return fault_sigs_; };

 protected:
  virtual void do_fault_injection();

  /** Creates the values that will be used in case of a fault
   *  Side effect: populates state2faultval_ and faultval2state_
   */
  virtual void create_fault_vals();

  TransitionSystem & fts_;
  FunctionalTransitionSystem faulty_fts_;

  smt::UnorderedTermMap state2faultsel_;
  smt::UnorderedTermMap faultsel2state_;

  smt::UnorderedTermMap state2faultval_;
  smt::UnorderedTermMap faultval2state_;

  smt::UnorderedTermMap state2faultsig_;
  smt::UnorderedTermMap faultsig2state_;

  smt::TermVec fault_sigs_;

  smt::UnorderedTermSet statevars_to_ignore_; ///< will not inject faults for these statevars
};

}  // namespace pono
