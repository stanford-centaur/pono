/*********************                                                        */
/*! \file nondet_fault_injector.cpp
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

#include "modifiers/fault_injector.h"

namespace pono {
class NonDetFaultInjector : public FaultInjector
{
 public:
  NonDetFaultInjector(FunctionalTransitionSystem & fts,
                      smt::UnorderedTermSet to_ignore = {})
      : FaultInjector(fts, to_ignore)
  {
    if (!fts.is_functional()) {
      throw PonoException("Expecting a functional transition system.");
    }
    do_fault_injection();
  }

 protected:
  void create_fault_vals() override;
};

}  // namespace pono
