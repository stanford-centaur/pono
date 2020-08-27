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
  Prover(const Property & p, smt::SolverEnum se);
  Prover(const Property & p, const smt::SmtSolver & s);
  Prover(const PonoOptions & opt, const Property & p, smt::SolverEnum se);
  Prover(const PonoOptions & opt, const Property & p, const smt::SmtSolver & s);

  virtual ~Prover();

  virtual void initialize();

  virtual ProverResult check_until(int k) = 0;

  bool witness(std::vector<smt::UnorderedTermMap> & out);

  ProverResult prove();
  
 protected:
  smt::SmtSolver solver_;
  smt::TermTranslator to_prover_solver_;
  const Property property_;
  const TransitionSystem &
      ts_;  ///< convenient reference to transition system in property
  const TransitionSystem &
      orig_ts_;  ///< reference to original TS before copied to new solver

  Unroller unroller_;

  int reached_k_;

  smt::Term bad_;

  PonoOptions options_;


  //TODO: can make some of COI functions/members private
  void print_bad_property_coi();

  void print_term_dfs(const smt::Term & term);

  void compute_term_coi(const smt::Term & term);

  void compute_coi();
  
  /* For static cone-of-influence analysis: 
     TermSets containing those state and input variables that appear
     in the term 'bad_' that represents the bad-state property. This
     information is used to rebuild the relevant parts of the
     transition system 'ts_' of the property. The (potentially)
     reduced and rebuilt transition system is stored in 'coi_ts_'. */
  smt::UnorderedTermSet statevars_in_coi_;
  smt::UnorderedTermSet inputvars_in_coi_;
  //  smt::UnorderedTermSet next_statevars_;
  TransitionSystem coi_ts_;
};
}  // namespace pono
