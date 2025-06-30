/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan, Florian Lonsing, Áron Ricardo Perez-Lopez
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief Base classes for algorithms used to prove LTL properties.
 **/

#pragma once

#include "core/prop.h"
#include "core/proverresult.h"
#include "core/ts.h"
#include "core/unroller.h"
#include "options/options.h"
#include "smt-switch/smt.h"

namespace pono {

class BaseProver
{
 public:
  virtual void initialize();

  virtual ProverResult prove();

  virtual ProverResult check_until(int k) = 0;

  virtual ~BaseProver() = default;

 protected:
  BaseProver(const TransitionSystem & ts, const smt::SmtSolver & solver);

  BaseProver(const TransitionSystem & ts,
             const smt::SmtSolver & solver,
             PonoOptions opt);

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

  /** Returns the reference of the interface ts, which is a copy of orig_ts but
   *  built using solver_. By default, the method returns a reference to ts_.
   *  The derived classes may be based on abstraction-refinement methods (e.g.
   *  IC3IA). In that case, the method would return the concrete ts.
   */
  virtual TransitionSystem & prover_interface_ts();

  smt::SmtSolver solver_;
  smt::TermTranslator to_prover_solver_;
  TransitionSystem orig_ts_;  ///< original TS before transferring to new solver
  TransitionSystem ts_;
  PonoOptions options_;
  Engine engine_ = Engine::NONE;

  bool initialized_ = false;
  int reached_k_ = -1;  ///< the last bound reached with no counterexamples

};  // class BaseProver

class SafetyProver : public BaseProver
{
 public:
  SafetyProver(const SafetyProperty & p,
               const TransitionSystem & ts,
               const smt::SmtSolver & s,
               PonoOptions opt = PonoOptions());

  void initialize() override;

  virtual bool witness(std::vector<smt::UnorderedTermMap> & out);

  /** Returns length of the witness
   *  this can be cheaper than actually computing the witness
   *  by default returns reached_k_+1, because reached_k_ was the
   *  last step that completed without finding a bug
   *  but some algorithms such as IC3 might need to follow the trace
   */
  virtual size_t witness_length() const;

  /** Gives a term representing an inductive invariant over current state
   * variables. Only valid if the property has been proven true. Only supported
   * by some engines
   */
  virtual smt::Term invar();

 protected:
  /** Default implementation for computing a witness
   *  Assumes that this engine is unrolling-based and that the solver
   *   state is currently satisfiable with a counterexample trace
   *  populates witness_
   *  @return true on success
   */
  virtual bool compute_witness();

  SafetyProperty orig_property_;  ///< original property before transferring

  Unroller unroller_;

  smt::Term bad_;

  // NOTE: both witness_ and invar_ use terms from the engine's solver

  std::vector<smt::UnorderedTermMap>
      witness_;  ///< populated by a witness if a CEX is found

  smt::Term invar_;  ///< populated with an invariant if the engine supports it
};  // class SafetyProver

using Prover [[deprecated("Use SafetyProver.")]] = SafetyProver;
}  // namespace pono
