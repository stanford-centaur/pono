/*********************                                                        */
/*! \file single_bit_fault_injector.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Modifies a transition system such that valid traces include
**        ones where a single bit can be flipped
**
**/

#pragma once

#include "modifiers/fault_injector.h"

namespace pono {
class SingleBitFaultInjector : public FaultInjector
{
 public:
  SingleBitFaultInjector(FunctionalTransitionSystem & fts,
                         smt::UnorderedTermSet to_ignore = {})
      : FaultInjector(fts, to_ignore)
  {
    if (!fts.is_functional()) {
      throw PonoException("Expecting a functional transition system.");
    }
    // TODO fix virtual calls in constructor!
    // can't call a virtual function in base class constructor
    // it will use the one from the base class, not the derived class
    // this is an issue since we use that here
    // maybe we should have an abstract base class and implement everything here
    do_fault_injection();
    constrain_to_single_fault();
  }

 protected:
  void create_fault_vals() override;

  // add constraints so that a fault can only occur at one time-step
  void constrain_to_single_fault();
};
}  // namespace pono
