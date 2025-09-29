/*********************                                                  */
/*! \file prop_monitor.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Adds a monitor for the property term and also replaces the property
**
**
**/

#pragma once

#include "core/ts.h"

namespace pono {

/** Adds a monitor state variable to the transition system for the given
 *  property. The monitor variable is initialized to true and holds the value of
 * `prop` at the previous time step. The transition system is modified in place.
 *
 *  The function can only be used when
 *  (1) `prop` contains a next-state variable (and `ts` is relational) or
 *  (2) the transition relation `ts.trans()` is right-total.
 *  This is to ensure that `ts` allows at least one further transition when
 * `prop` is falsified (recall that the monitor has a one-step delay).
 *
 *  NOTE: only works for safety properties currently.
 *
 *  @param ts the transition system to modify
 *  @param prop the property term to monitor
 *  @return the monitor state variable term
 */
smt::Term add_prop_monitor(TransitionSystem & ts, const smt::Term & prop);

}  // namespace pono
