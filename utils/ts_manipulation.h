/*********************                                                        */
/*! \file ts_manipulation.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Useful functions for creating and manipulating transition systems
**
**
**/
#pragma once

#include "core/ts.h"

namespace pono {

// returns an empty system over the given solver
TransitionSystem create_fresh_ts(bool functional,
                                 const smt::SmtSolver & solver);

}  // namespace pono
