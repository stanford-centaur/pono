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

/** returns true iff l is a literal
 *  e.g. either a boolean symbolic constant or its negation
 *  NOTE will return false for nested negations, i.e. (not (not (not l)))
 *  @param l the term to check
 *  @param boolsort a boolean sort from the corresponding solver
 *         this way sort aliasing solvers are still supported
 *  @return true iff l is a literal
 */
bool is_lit(const smt::Term & l, const smt::Sort & boolsort);

/** returns all the free (not bound to a quantifier) symbols
 *  in a Term
 *  Note: includes uninterpreted functions
 */
smt::UnorderedTermSet get_free_symbols(const smt::Term & term);

/** Extract all predicates from a term
 *  Traverses all the subterms of term and adds any predicates to out
 *  @param the solver to use (needed for building new terms when processing
 * ITEs)
 *  @param term the term to traverse
 *  @param out the set to add to
 *  @param include_symbols if set to true, will include boolean symbols
 *         (0-arity predicates), otherwise will ignore them
 *
 *  TEMPORARY RESTRICTION
 *    doesn't work with boolector as backend solver if there are uninterpreted
 * functions the smt-switch btor backend doesn't support getting the sort from a
 * UF
 */
void get_predicates(const smt::SmtSolver & solver,
                    const smt::Term & term,
                    smt::UnorderedTermSet & out,
                    bool include_symbols = false);

}  // namespace pono
