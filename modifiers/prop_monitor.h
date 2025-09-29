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

#include "core/rts.h"
#include "core/ts.h"
#include "utils/logger.h"

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
smt::Term add_prop_monitor(TransitionSystem & ts, const smt::Term & prop)
{
  // only in debug mode as right-total check can be expensive
  assert(!ts.no_next(prop) || ts.is_right_total());
  logger.log(1, "Adding a monitor for the property");

  // checks for functional TS
  if (ts.is_functional() && !ts.no_next(prop)) {
    throw PonoException(
        "Cannot use next in property of a functional transition system.");
  }

  smt::Term monitor;
  size_t id = 0;
  while (true) {
    try {
      monitor = ts.make_statevar("_monitor_" + std::to_string(id),
                                 ts.make_sort(smt::BOOL));
      break;
    }
    catch (SmtException & e) {
      ++id;
    }
  }
  assert(monitor);

  // monitor starts true
  ts.constrain_init(monitor);

  if (ts.no_next(prop)) {
    ts.assign_next(monitor, prop);
  } else {
    assert(!ts.is_functional());
    RelationalTransitionSystem & rts =
        static_cast<RelationalTransitionSystem &>(ts);
    rts.constrain_trans(rts.make_term(smt::Equal, rts.next(monitor), prop));
  }

  return monitor;
}

}  // namespace pono
