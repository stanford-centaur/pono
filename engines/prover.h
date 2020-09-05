/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan, Florian Lonsing
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#pragma once

#include "core/prop.h"
#include "core/proverresult.h"
#include "core/ts.h"
#include "core/unroller.h"
#include "options/options.h"

#include "smt-switch/smt.h"

namespace pono {
class Prover
{
 public:
  Prover(Property & p, smt::SolverEnum se);
  Prover(Property & p, const smt::SmtSolver & s);
  Prover(const PonoOptions & opt, Property & p, smt::SolverEnum se);
  Prover(const PonoOptions & opt, Property & p, const smt::SmtSolver & s);

  virtual ~Prover();

  virtual void initialize();

  virtual ProverResult prove();

  virtual ProverResult check_until(int k) = 0;

  virtual bool witness(std::vector<smt::UnorderedTermMap> & out);

  /** Gives a term representing an inductive invariant over current state
   * variables. Only valid if the property has been proven true. Only supported
   * by some engines
   */
  virtual smt::Term invar();

 protected:
  /** Take a term from the Prover's solver
   *  to the original transition system's solver
   *  as a particular SortKind
   *  Note: they could be the same solver but aren't necessarily
   *  @param t a term from the Prover's solver
   *  @param sk the SortKind to cast with
   *  @return a term in the original TS's solver
   */
  smt::Term to_orig_ts(smt::Term t, smt::SortKind sk);

  /** Take a term from the Prover's solver
   *  to the original transition system's solver
   *  Note: they could be the same solver but aren't necessarily
   *  @param t a term from the Prover's solver
   *  @return a term in the original TS's solver
   */
  smt::Term to_orig_ts(smt::Term t);

  smt::SmtSolver solver_;
  smt::TermTranslator to_prover_solver_;
  Property property_;
  TransitionSystem &
      ts_;  ///< convenient reference to transition system in property
  TransitionSystem &
      orig_ts_;  ///< reference to original TS before copied to new solver

  Unroller unroller_;

  int reached_k_;

  smt::Term bad_;

  PonoOptions options_;

 private:
  /* Cone-of-influence analysis. */

  /* Debugging helper functions. */
  void print_coi_info();
  void print_term_dfs(const smt::Term & term);
  /* Key functions. */
  void compute_coi();
  void collect_coi_term(smt::UnorderedTermSet & set, const smt::Term & term);
  void compute_coi_trans_constraints();
  void compute_term_coi(const smt::Term & term,
                        smt::UnorderedTermSet & new_coi_state_vars,
                        smt::UnorderedTermSet & new_coi_input_vars);
  void compute_coi_next_state_funcs();

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

  /* End: cone-of-influence analysis. */
};
}  // namespace pono
