/*********************                                                  */
/*! \file prophecy_modifier.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Class for adding prophecy variables that predict a value some
**        constant number of steps before a property violation.
**/
#pragma once

#include <utility>

#include "core/ts.h"
#include "modifiers/history_modifier.h"

namespace pono {

class ProphecyModifier
{
 public:
  ProphecyModifier(TransitionSystem & ts);

  /** Returns a prophecy variable predicting the target delay steps
   *  before a property violation (if the property can be violated)
   *  will create new variables and update the transition system as needed
   *  @param target a current state variable to target for a prophecy variable
   *  @param delay the amount of delay to introduce for the prophecy
   *  @param prop the current property
   *  @return the prophecy variable and the updated target
   *          (a history variable if delay is non-zero)
   *          this should be used to update the property. i.e.
   *          if P was the original property, it should now be
   *          proph=target -> P
   */
  std::pair<smt::Term, smt::Term> get_proph(const smt::Term & target,
                                            size_t delay);

 protected:
  TransitionSystem & ts_;
  const smt::SmtSolver solver_;
  HistoryModifier hm_;

  // maps current state variables to a list of prophecy variables
  // where the index corresponds to the delay of the target (-1)
  std::unordered_map<smt::Term, smt::TermVec> proph_vars_;
};

}  // namespace pono
