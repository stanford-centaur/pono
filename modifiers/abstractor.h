/*********************                                                  */
/*! \file abstractor.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Generic base class for abstractors -- classes that can
**        create an abstract version of a transition system
**
**/
#pragma once

#include "core/fts.h"
#include "core/rts.h"

namespace pono {

class Abstractor
{
 public:
  Abstractor(const TransitionSystem & ts)
      : conc_ts_(ts), abs_ts_(conc_ts_.solver())
  {
    // TODO: handle functional / relational correctly
    //       issue is constructor and copy assignment not working well
    //       the abs_ts_ will be initialized with a default solver
    //       resulted in a segfault in CVC4
    //       currently just always using a relational system
  }

  virtual ~Abstractor() {}

  /** Returns the abstraction of a concrete term
   *  @param the concrete term to abstract
   *  @return the abstracted term
   *  This is an empty implementation. Derived classes will implement this.
   */
  virtual smt::Term abstract(smt::Term & t)
  {
    throw PonoException("Abstractor base class does not implement methods.");
  }

  /** Returns the concretization of an abstract term
   *  @param the abstract term to concretize
   *  @return the concrete version of the term
   *  This is a NOP implementation. Derived classes will implement this.
   */
  virtual smt::Term concrete(smt::Term & t)
  {
    throw PonoException("Abstractor base class does not implement methods.");
  }

  /** Getter for the abstracted transition system
   *  This is intentionally a non-const reference because often
   *  the abstraction will be refined. This way, it does not
   *  need to be copied before being refined.
   *  @return a reference to the abstracted system
   */
  TransitionSystem & abs_ts() { return abs_ts_; }

 protected:
  /** Perform the abstraction
   *  This should populate the abstraction and concretization caches
   *  This is a NOP implementation. Derived classes will implement this.
   *  Should be called in derived class constructor.
   */
  virtual void do_abstraction(){};

  const TransitionSystem & conc_ts_;
  RelationalTransitionSystem abs_ts_;

  smt::UnorderedTermMap abstraction_cache_;
  smt::UnorderedTermMap concretization_cache_;
};

}  // namespace pono
