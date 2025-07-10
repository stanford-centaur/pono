/*********************                                                        */
/*! \file cegar_ops_uf.h
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief A simple CEGAR loop that abstracts operators with uninterpreted
**        functions and refines by concretizing the operators
**
**/

#pragma once

#include "core/unroller.h"
#include "engines/cegar.h"
#include "modifiers/ops_abstractor.h"

namespace pono {

template <class Prover_T>
class CegarOpsUf : public CEGAR<Prover_T>
{
  typedef CEGAR<Prover_T> super;

 public:
  CegarOpsUf(const SafetyProperty & p,
             const TransitionSystem & ts,
             const smt::SmtSolver & solver,
             PonoOptions opt = PonoOptions());

  void set_ops_to_abstract(const smt::UnorderedOpSet & ops_to_abstract);

  void set_min_bitwidth(unsigned long w) { oa_.set_min_bitwidth(w); };

  void initialize() override;

  ProverResult check_until(int k) override;

 protected:
  void cegar_abstract() override;

  bool cegar_refine() override;

  void refine_subprover_ts_base(const smt::UnorderedTermSet & axioms,
                                bool skip_init);
  void refine_subprover_ts(const smt::UnorderedTermSet & axioms,
                           bool skip_init);

  TransitionSystem & prover_interface_ts() override { return conc_ts_; }

  TransitionSystem conc_ts_;
  TransitionSystem & prover_ts_;

  OpsAbstractor oa_;

  // solver and associated infrastructure for
  // unrolling based refinement
  smt::SmtSolver cegopsuf_solver_;
  smt::TermTranslator to_cegopsuf_solver_;
  smt::TermTranslator from_cegopsuf_solver_;
  TransitionSystem cegopsuf_ts_;
  Unroller cegopsuf_un_;
  smt::Term cegopsuf_bad_;

  smt::UnorderedTermMap cegopsuf_labels_;  // labels for each abstract uf
};

}  // namespace pono
