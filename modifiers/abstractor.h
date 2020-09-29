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
  Abstractor(const TransitionSystem & conc_ts, TransitionSystem & abs_ts)
      : conc_ts_(conc_ts), abs_ts_(abs_ts)
  {
    if (abs_ts_.is_functional() && !conc_ts_.is_functional()) {
      // this might not be true in every case
      // but playing it safe now
      // this is at least true for all our abstraction ideas so far
      throw PonoException(
          "Cannot abstract a relational system with a functional system");
    }
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

  // getters
  const TransitionSystem & conc_ts() const { return conc_ts_; };

  // TODO: consider making this const and just copying the system
  /** The abstract transition system is non-const
   *  so that it can be refined in-place.
   */
  TransitionSystem & abs_ts() const { return abs_ts_; };

 protected:
  /** Perform the abstraction
   *  This should populate the abstraction and concretization caches
   *  This is a NOP implementation. Derived classes will implement this.
   *  Should be called in derived class constructor.
   */
  virtual void do_abstraction(){};

  const TransitionSystem & conc_ts_;
  TransitionSystem & abs_ts_;

  smt::UnorderedTermMap abstraction_cache_;
  smt::UnorderedTermMap concretization_cache_;
};

}  // namespace pono
