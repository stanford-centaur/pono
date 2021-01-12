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

/** returns true iff t is a predicate
 *  @param t the term to check
 *  @param boolsort a boolean sort from the corresponding solver
 *         this way sort aliasing solvers are still supported
 *  @return true iff t is a predicate
 */
bool is_predicate(const smt::Term & t, const smt::Sort & boolsort);

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
 *  @param search_subterms if true continues to look for predicates
 *         in children of a predicate. This can happen when there are ITEs
 *  @param split_ites if true then don't include ITEs in predicates, always
 *         consider the different cases instead
 *
 *  TEMPORARY RESTRICTION
 *    doesn't work with boolector as backend solver if there are uninterpreted
 * functions the smt-switch btor backend doesn't support getting the sort from a
 * UF
 */
void get_predicates(const smt::SmtSolver & solver,
                    const smt::Term & term,
                    smt::UnorderedTermSet & out,
                    bool include_symbols = false,
                    bool search_subterms = false,
                    bool split_ites = false);

}  // namespace pono
