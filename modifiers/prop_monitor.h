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
#include "smt-switch/smt.h"

namespace pono {

smt::Term add_prop_monitor(TransitionSystem & ts, const smt::Term & prop);

}  // namespace pono
