/*********************                                                        */
/*! \file term_analysis.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Useful functions for term analysis.
**
**
**/

#pragma once

#include "smt-switch/smt.h"

namespace pono {

/** adds all the free (not bound to a quantifier) symbols
 *  in a Term to an existing set
 *  Note: includes uninterpreted functions
 */
void get_free_symbols(const smt::Term & term,
                      smt::UnorderedTermSet & out_symbols);

/** returns all the free (not bound to a quantifier) symbols
 *  in a Term
 *  Note: includes uninterpreted functions
 */
smt::UnorderedTermSet get_free_symbols(const smt::Term & term);
}  // namespace pono
