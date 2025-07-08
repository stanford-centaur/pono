/*********************                                                        */
/*! \file control_signals.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Useful functions for adding control signal semantics to a system.
**        For example, toggling a clock or setting a reset sequence.
**
**/
#pragma once

#include "core/ts.h"

namespace pono {

/** Toggle a clock symbol. Currently only supports starting at 0 and
 *  toggling every step. This can be extended in the future.
 *
 *  @param ts the TransitionSystem to modify
 *  @param clock_symbol the clock to toggle
 *         assumed to be a boolean or bit-vector of size one input or current
 * var where posedge is when it's evaluated to true/bv1
 */
void toggle_clock(TransitionSystem & ts, const smt::Term & clock_symbol);

/** Holds a reset signal active for reset_bnd steps starting in the first state.
 *  Returns the condition to guard a property with to not check
 *  it until after the reset sequence has ended.
 *
 *  @param ts the transition system to modify
 *  @param reset_symbol the reset signal (assumed to be in the transition
 * system) assumed to be a boolean or bit-vector of size one input or current
 * var where reset is active when it's evaluated to true/bv1
 *  @param reset_bnd how many steps to hold reset active
 *  @return a term which is true in ts when the reset sequence is done
 */
smt::Term add_reset_seq(TransitionSystem & ts,
                        const smt::Term & reset_symbol,
                        unsigned long reset_bnd = 1);

}  // namespace pono
