/*********************                                                  */
/*! \file coi.h
** \verbatim
** Top contributors (to current version):
**   Florian Lonsing, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Class for performing cone of influence reduction
**
**
**/

#pragma once

#include "core/ts.h"

namespace pono {
class ConeOfInfluence
{
 public:
  /** This class computes the cone of influence on construction
   *  @param ts the transition system to modify
   *  @param to_keep terms in the transition system that need to be kept
   *  The cone-of-influence will keep all the variables from terms in
   *    to_keep and any variables that influence those variables.
   */
  ConeOfInfluence(TransitionSystem & ts,
                  const smt::TermVec & to_keep,
                  const smt::TermVec & to_remove = {},
                  int verbosity = 1);

 protected:
  /* Debugging helper functions. */
  void print_coi_info(const smt::TermVec & to_keep);
  void print_term_dfs(const smt::Term & term);
  /* Key functions. */
  void compute_coi(const smt::TermVec & to_keep);
  void collect_coi_term(smt::UnorderedTermSet & set, const smt::Term & term);
  void compute_coi_trans_constraints();
  void compute_term_coi(const smt::Term & term,
                        smt::UnorderedTermSet & new_coi_state_vars,
                        smt::UnorderedTermSet & new_coi_input_vars);
  void compute_coi_next_state_funcs();

  TransitionSystem & ts_;
  int verbosity_;
  smt::UnorderedTermSet to_remove_;
  /* TermSets containing those state and input variables that appear
     in the term 'bad_' that represents the bad-state property. This
     information is used to rebuild the transition relation of the
     transition system 'ts_' of the property. */
  smt::UnorderedTermSet statevars_in_coi_;
  smt::UnorderedTermSet inputvars_in_coi_;
  /* Set of terms already visited in COI analysis. */
  smt::UnorderedTermSet coi_visited_terms_;
  unsigned int orig_num_statevars_;
  unsigned int orig_num_inputvars_;
};
}  // namespace pono
