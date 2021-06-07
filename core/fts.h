/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#pragma once

#include <unordered_map>

#include "core/ts.h"
#include "utils/exceptions.h"

namespace pono {

class FunctionalTransitionSystem : public TransitionSystem
{
  typedef TransitionSystem super;

 public:
  FunctionalTransitionSystem() : TransitionSystem() { functional_ = true; }

  FunctionalTransitionSystem(const smt::SmtSolver & s) : TransitionSystem(s)
  {
    functional_ = true;
  }

  FunctionalTransitionSystem(const TransitionSystem & other_ts,
                             smt::TermTranslator & tt)
      : TransitionSystem(other_ts, tt)
  {
    functional_ = true;
  }

  FunctionalTransitionSystem(const TransitionSystem & other_ts)
      : TransitionSystem(other_ts)
  {
    if (!functional_) {
      throw PonoException("Copied relational system to functional one");
    }
  }

 protected:

  // helpers and checkers

  // overloaded
  bool known_symbols(const smt::Term & term) const override;

  // Note: this method is not exposed, because the user should not be able to
  //       build terms with next states
  /* Replace all current states by their next-state updates, functionally */
  smt::Term to_next_func(const smt::Term & term);
};

}  // namespace pono
