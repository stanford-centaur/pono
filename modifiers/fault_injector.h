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
  FaultInjector(FunctionalTransitionSystem & fts,
                smt::UnorderedTermSet to_ignore = {})
      : fts_(fts),
        faulty_fts_(fts),
        statevars_to_ignore_(to_ignore),
        initialized_(false)
  {
    if (!fts.is_functional()) {
      throw PonoException("Can only do fault injection on functional systems.");
    }
  };

  virtual ~FaultInjector(){};

  /** Inject faults in the transition system if that hasn't already been done,
      and return the resulting system
      @return a TransitionSystem with injected faults.
   */
  FunctionalTransitionSystem & faulty_transition_system()
  {
    if (!initialized_) {
      do_fault_injection();
    }
    return faulty_fts_;
  };

  bool is_initialized() const { return initialized_; };

  smt::TermVec fault_sigs() const
  {
    if (!initialized_) {
      throw PonoException(
          "Cannot get fault signals without calling faulty_transition_system "
          "first");
    }
    return fault_sigs_;
  };

 protected:
  virtual void do_fault_injection();

  /** Creates the values that will be used in case of a fault
   *  Side effect: populates state2faultval_ and faultval2state_
   */
  virtual void create_fault_vals() = 0;

  FunctionalTransitionSystem & fts_;
  FunctionalTransitionSystem faulty_fts_;

  smt::UnorderedTermSet
      statevars_to_ignore_;  ///< will not inject faults for these statevars

  bool initialized_;  ///< set to true once faults have been injected

  smt::UnorderedTermMap state2faultsel_;
  smt::UnorderedTermMap faultsel2state_;

  smt::UnorderedTermMap state2faultval_;
  smt::UnorderedTermMap faultval2state_;

  smt::UnorderedTermMap state2faultsig_;
  smt::UnorderedTermMap faultsig2state_;

  smt::TermVec fault_sigs_;

};

}  // namespace pono
