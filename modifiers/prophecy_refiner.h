/*********************                                                  */
/*! \file prophecy_refiner.h
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
#include "modifiers/history_refiner.h"

namespace pono {

class ProphecyRefiner
{
 public:
  ProphecyRefiner(TransitionSystem & ts);

  /** Returns a prophecy variable predicting the target delay steps
   *  before a property violation (if the property can be violated)
   *  will create new variables and update the transition system as needed
   *  @param target a current state variable to target for a prophecy variable
   *  @param delay the amount of delay to introduce for the prophecy
   *  @param prop the current property
   *  @return the prophecy variable and the updated property
   */
  std::pair<smt::Term, smt::Term> get_proph(const smt::Term & target,
                                            size_t delay,
                                            const smt::Term & prop);

 protected:
  TransitionSystem & ts_;
  const smt::SmtSolver solver_;
  HistoryRefiner hr_;

  // maps current state variables to a list of prophecy variables
  // where the index corresponds to the delay of the target (-1)
  std::unordered_map<smt::Term, smt::TermVec> proph_vars_;
};

}  // namespace pono
