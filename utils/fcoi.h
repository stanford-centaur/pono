/*********************                                                  */
/*! \file static_coi.h
** \verbatim
** Top contributors (to current version):
**   Florian Lonsing, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Class for finding symbols in the cone-of-influence of terms
**        It traverses next-state functions and constraints to determine
**        all symbols that can influence those terms
**
**/

#pragma once

#include "core/ts.h"
#include "utils/logger.h"

namespace pono {
class FunctionalConeOfInfluence
{
 public:
  /** This class computes the cone of influence on construction
   *  @param ts the transition system to modify
   *  @param to_keep terms in the transition system that need to be kept
   *  The cone-of-influence will keep all the variables from terms in
   *    to_keep and any variables that influence those variables.
   */
  FunctionalConeOfInfluence(const TransitionSystem & ts, int verbosity = 1);

  // TODO if this is used incrementally, consider caching
  //      some info. could be optional

  /** Compute the cone of influence for terms
   *  @param terms - a vector of important terms which
   *  after running, can access that statevars and input
   *  vars which fall into the coi with
   *     statevars_in_coi() and
   *     inputvars_in_coi(), respectively
   *
   *  those are only available until the next call
   *  to compute_coi, which automatically clears
   *  those data structures before re-running
   */
  void compute_coi(const smt::TermVec & terms);

  const smt::UnorderedTermSet & statevars_in_coi() const
  {
    return statevars_in_coi_;
  }

  const smt::UnorderedTermSet & inputvars_in_coi() const
  {
    return inputvars_in_coi_;
  }

 protected:
  /* Helper functions */

  /** Clears all internal data structures.
   *  Called automatically at beginning of compute_coi
   */
  void clear();

  /* Debugging helper functions. */
  void print_coi_info(const smt::TermVec & terms);
  void print_term_dfs(const smt::Term & term);

  /* Key functions. */
  void compute_coi_trans_constraints();
  void compute_term_coi(const smt::Term & term,
                        smt::UnorderedTermSet & new_coi_state_vars,
                        smt::UnorderedTermSet & new_coi_input_vars);
  void compute_coi_next_state_funcs();

  const TransitionSystem & ts_;
  int verbosity_;

  Log local_logger_;  ///< local instance of a logger to respect this verbosity

  /* TermSets containing those state and input variables that appear
     in the term 'bad_' that represents the bad-state property. This
     information is used to rebuild the transition relation of the
     transition system 'ts_' of the property. */
  smt::UnorderedTermSet statevars_in_coi_;
  smt::UnorderedTermSet inputvars_in_coi_;
  /* Set of terms already visited in COI analysis. */
  smt::UnorderedTermSet coi_visited_terms_;
};
}  // namespace pono
