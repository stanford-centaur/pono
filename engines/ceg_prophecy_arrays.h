/*********************                                                        */
/*! \file ceg_prophecy_arrays.h
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief An implementation of Counter-Example Guided Prophecy for array
**        model checking. It is parameterized by an underlying model checking
**        procedure which need not handle arrays (only UF). However, a common
**        instantiation is with an IC3-style procedure, in which case we
**        often refer to this algorithm as "prophic3".
**
**/

#pragma once

#include "core/unroller.h"
#include "engines/cegar.h"
#include "modifiers/array_abstractor.h"
#include "modifiers/prophecy_modifier.h"
#include "options/options.h"
#include "refiners/array_axiom_enumerator.h"

namespace pono {

template <class Prover_T>
class CegProphecyArrays : public CEGAR<Prover_T>
{
  typedef CEGAR<Prover_T> super;

 public:
  CegProphecyArrays(const SafetyProperty & p,
                    const TransitionSystem & ts,
                    const smt::SmtSolver & solver,
                    PonoOptions opt = PonoOptions());

  // HACK override prove so that it calls prove
  // in the underlying model checker instead of check_until
  // especially important if msat-ic3ia is the underlying engine
  // because it only supports prove
  ProverResult prove() override;

  ProverResult check_until(int k) override;

  void initialize() override;

  size_t witness_length() const override
  {
    // regardless of super class want this to be the witness length
    return reached_k_ + 1;
  }

 protected:
  TransitionSystem conc_ts_;
  TransitionSystem & abs_ts_;
  Unroller abs_unroller_;
  ArrayAbstractor aa_;
  ArrayAxiomEnumerator aae_;
  ProphecyModifier pm_;

  int reached_k_;  ///< local variable to check the length of BMC refinement run

  size_t num_added_axioms_;  ///< set by refine to the number of added axioms

  smt::UnorderedTermMap labels_;  ///< labels for unsat core minimization

  smt::UnorderedTermSet
      important_vars_;  ///< important variables
                        ///< useful for IC3IA to prioritize predicates

  TransitionSystem & prover_interface_ts() override { return conc_ts_; }

  void cegar_abstract() override;
  bool cegar_refine() override;

  // helpers
  smt::Term get_bmc_formula(size_t b);

  /** Unsat core based axiom reduction
   *  @param abs_bmc_formula the trace formula
   *  @param consec_ax consecutive axioms over transition system variables
   *         axioms are dropped from this set if they are not needed
   */
  void reduce_consecutive_axioms(const smt::Term & abs_bmc_formula,
                                 smt::UnorderedTermSet & consec_ax);

  /** Unsat core based axiom reduction
   *  @param abs_bmc_formula the trace formula
   *         Note: consecutive axioms need to be added to the trace
   *               formula outside of this function
   *  @param nonconsec_ax nonconsecutive axioms over unrolled variables
   *  @return a subset of the nonconsecutive axioms that are sufficient
   *    to rule out the trace formula
   */
  AxiomVec reduce_nonconsecutive_axioms(const smt::Term & abs_bmc_formula,
                                        const AxiomVec & nonconsec_ax);

  /** Lookup or create a label for t
   *  Uses and modifies labels_
   *  @param t a boolean term to create a label for
   *  @return the label
   */
  smt::Term label(const smt::Term & t);

  void refine_ts(const smt::UnorderedTermSet & consecutive_axioms);

  void refine_subprover_ts(const smt::UnorderedTermSet & consecutive_axioms);

  void add_important_var(const smt::Term & v);
};

}  // namespace pono
