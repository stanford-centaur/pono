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
  }

  typedef FaultInjector super;

 protected:
  // override functions that are called to inject faults
  void do_fault_injection() override;
  void create_fault_vals() override;

  // add constraints so that a fault can only occur at one time-step
  void constrain_to_single_fault();
};
}  // namespace pono
