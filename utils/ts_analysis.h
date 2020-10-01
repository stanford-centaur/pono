/*********************                                                        */
/*! \file ts_analysis.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Useful functions for analyzing transition systems.
**
**
**/
#pragma once

#include "smt-switch/smt.h"

#include "core/ts.h"

namespace pono {

/** Check if a given term over current state variables
 *  is an inductive invariant for the given system
 *  and property
 *  @param ts the transition system
 *  @param prop the term representing the property
 *  @param invar the term representing the invariant
 *  @return true iff invar is an inductive invariant that
 *          guarantees prop holds
 */
bool check_invar(const TransitionSystem & ts,
                 const smt::Term & prop,
                 const smt::Term & invar);

}  // namespace pono
