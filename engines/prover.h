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

/** enum for communicating result of a refinement step
 *  only used for algorithms that use abstraction refinement
 */
enum RefineResult
{
  REFINE_NONE = 0,  // no refinement necessary (e.g. concrete)
  REFINE_SUCCESS,   // refinement successful
  REFINE_FAIL       // failed to refine
};

class Prover
{
 public:
  Prover(const Property & p, const TransitionSystem & ts,
         const smt::SmtSolver & s,
         PonoOptions opt = PonoOptions());

  virtual ~Prover();

  virtual void initialize();

  virtual ProverResult prove();

  virtual ProverResult check_until(int k) = 0;

  virtual bool witness(std::vector<smt::UnorderedTermMap> & out);

  /** Gives a term representing an inductive invariant over current state
   * variables. Only valid if the property has been proven true. Only supported
   * by some engines
   */
  smt::Term invar();

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

  /** Default implementation for computing a witness
   *  Assumes that this engine is unrolling-based and that the solver
   *   state is currently satisfiable with a counterexample trace
   *  populates witness_
   *  @return true on success
   */
  bool compute_witness();

  bool initialized_;

  smt::SmtSolver solver_;
  smt::TermTranslator to_prover_solver_;

  Property orig_property_; ///< original property before copied to new solver
  TransitionSystem orig_ts_;  ///< original TS before copied to new solver

  TransitionSystem ts_;

  Unroller unroller_;

  int reached_k_;

  //TODO: add engine enum

  smt::Term bad_;

  PonoOptions options_;

  Engine engine_;

  // NOTE: both witness_ and invar_ are use terms from the engine's solver

  std::vector<smt::UnorderedTermMap> witness_; ///< populated by a witness if a CEX is found

  smt::Term invar_; ///< populated with an invariant if the engine supports it

};
}  // namespace pono
