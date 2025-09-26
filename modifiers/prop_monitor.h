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
 *  property. The monitor variable is true iff the property has held at all
 *  previous time steps. The transition system is modified in place.
 *
 *  NOTE: only works for safety properties currently.
 *
 *  NOTE: cannot be used on functional transition systems with constraints (more
 *    precisely, non-right-total transition relation).
 *
 *  @param ts the transition system to modify
 *  @param prop the property term to monitor
 *  @return the monitor state variable term
 */
smt::Term add_prop_monitor(TransitionSystem & ts, const smt::Term & prop)
{
  logger.log(1, "Adding a monitor for the property");

  // checks for functional TS
  if (ts.is_functional()) {
    if (!ts.no_next(prop)) {
      throw PonoException(
          "Cannot use next in property of a functional transition system.");
    }
    if (!ts.constraints().empty()) {
      // TODO: more precise right-total check
      // (https://github.com/stanford-centaur/pono/pull/469)
      throw PonoException(
          "Cannot add monitor on functional transition systems "
          "with constraints.");
    }
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

  smt::SmtSolver s = ts.solver();
  // monitor is set to false if prop is ever false:
  // monitor' = prop && monitor
  smt::Term monitor_next = s->make_term(smt::And, prop, monitor);
  if (ts.is_functional()) {
    assert(ts.no_next(prop));
    ts.assign_next(monitor, monitor_next);
  } else {  // relational TS
    RelationalTransitionSystem & rts =
        static_cast<RelationalTransitionSystem &>(ts);
    if (rts.no_next(prop)) {
      rts.assign_next(monitor, monitor_next);
    } else {
      rts.constrain_trans(
          rts.make_term(smt::Equal, rts.next(monitor), monitor_next));
    }
    // ensure that if prop is false, there will always be a next state
    rts.set_trans(
        s->make_term(smt::Or, s->make_term(smt::Not, prop), rts.trans()));
  }

  return monitor;
}

}  // namespace pono
