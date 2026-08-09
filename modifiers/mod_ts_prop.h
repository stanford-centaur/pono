/*********************                                                  */
/*! \file mod_ts_prop.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Functions for modifying ts and property
**
**
**/

#pragma once

#include "core/ts.h"

namespace pono {

/** Creates new state variables to represent init and trans
 *  and adds an extra transition from the pseudo-initial state
 *  and to the property monitor
 *
 *  NOTE: doing this requires making the ts relational
 *
 *  @param ts the transition system to modify
 *  @param prop the property to modify
 *  @return the updated transition system
 *          the returned system is relational even if ts is functional
 *  Updates the prop in-place
 */
TransitionSystem pseudo_init_and_prop(TransitionSystem & ts, smt::Term & prop);

// optimization to assume the property in the pre-state
// although confusing, this is sound as long as you always check
// for a property violation over the next-state variables
void prop_in_trans(TransitionSystem & ts, const smt::Term & prop);

/** Returns a new transition system with the given input variables
 *  promoted to state variables with no update (e.g. implicit inputs).
 *  This pass automatically re-evaluates constraints which might need to be
 *  added differently over state variables vs input variables
 *  @param ts the transition system to promote input variables in
 *  @param ivs_to_promote input variables to promote
 *  @return a transition system that is semantically the same but
 *          with no input variables (only implicit inputs)
 */
TransitionSystem promote_inputvars(
    const TransitionSystem & ts, const smt::UnorderedTermSet & ivs_to_promote);

/** Returns a new transition system with *all* input variables in the transition
 *  system promoted to state variables with no update (e.g. implicit inputs).
 *  This pass automatically re-evaluates constraints which might need to be
 *  added differently over state variables vs input variables
 *  @param ts the transition system to promote input variables in
 *  @return a transition system that is semantically the same but
 *          with no input variables (only implicit inputs)
 */
TransitionSystem promote_inputvars(const TransitionSystem & ts);

/**
 * @brief Instrument the transition system to make it right-total.
 *
 * An auxiliary state variable `valid` is added to the system such that
 * - `valid` is true iff the original transition relation holds
 *   in all previous time steps
 * - the new property becomes `valid => prop`
 *
 * A new transition system and a new property term will be
 * created and replace `ts` and `prop`, respectively
 *
 * @param ts the transition system to modify
 * @param prop the safety property
 */
void make_trans_total(TransitionSystem & ts, smt::Term & prop);

}  // namespace pono
