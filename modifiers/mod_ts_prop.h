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
void prop_in_trans(TransitionSystem & ts, const smt::Term & prop)
{
  // NOTE: CRUCIAL that we pass false here
  // cannot add to init or the next states
  // passing false prevents that
  ts.add_constraint(prop, false);
}

}  // namespace pono
