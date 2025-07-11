/*********************                                                        */
/*! \file make_provers.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Utility functions for creating provers.
**
**
**/

#pragma once

#include "engines/prover.h"

namespace pono {

std::vector<Engine> all_engines();

std::shared_ptr<SafetyProver> make_prover(Engine e,
                                          const SafetyProperty & p,
                                          const TransitionSystem & ts,
                                          const smt::SmtSolver & slv,
                                          PonoOptions opts = PonoOptions());

std::shared_ptr<SafetyProver> make_ceg_proph_prover(
    Engine e,
    const SafetyProperty & p,
    const TransitionSystem & ts,
    const smt::SmtSolver & slv,
    PonoOptions opts = PonoOptions());

std::shared_ptr<SafetyProver> make_cegar_values_prover(
    Engine e,
    const SafetyProperty & p,
    const TransitionSystem & ts,
    const smt::SmtSolver & slv,
    PonoOptions opts = PonoOptions());

std::shared_ptr<SafetyProver> make_cegar_bv_arith_prover(
    Engine e,
    const SafetyProperty & p,
    const TransitionSystem & ts,
    const smt::SmtSolver & slv,
    PonoOptions opts = PonoOptions(),
    const smt::UnorderedOpSet & ops_to_abstract = { { smt::BVMul,
                                                      smt::BVUdiv,
                                                      smt::BVSdiv,
                                                      smt::BVUrem,
                                                      smt::BVSrem,
                                                      smt::BVSmod } });

}  // namespace pono
