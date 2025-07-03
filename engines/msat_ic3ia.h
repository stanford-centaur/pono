/*********************                                                        */
/*! \file msat_ic3ia.h
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief A backend using the open-source IC3IA implementation from
**        Alberto Griggio. Available here:
**        https://es-static.fbk.eu/people/griggio/ic3ia/index.html
**
**/

#pragma once

#include "engines/prover.h"
#include "ic3ia/ic3.h"

namespace pono {

class MsatIC3IA : public SafetyProver
{
 public:
  MsatIC3IA(const SafetyProperty & p,
            const TransitionSystem & ts,
            const smt::SmtSolver & solver,
            PonoOptions opt = PonoOptions());

  typedef SafetyProver super;

  ProverResult prove() override;

  // TODO: might have to ignore this and only support prove
  // ic3ia doesn't have a way to only check up to some number of frames
  // might need to make prove virtual as well
  ProverResult check_until(int k) override;

 protected:
  /** Recover a witness from an IC3IA object
   *  populates the witness_ member variable with terms from this engine's
   * solver
   *  @param env the msat_env used by the IC3IA backend
   *  @param ic3 the backend ic3ia engine that was used -- should have returned
   * false
   *  @param to_ts_solver a TermTranslator that moves terms from the fresh
   * solver for ic3ia to terms of solver_
   *  @return true on success
   */
  bool compute_witness(msat_env env,
                       ic3ia::IC3 & ic3,
                       smt::TermTranslator & to_ts_solver);
};

}  // namespace pono
