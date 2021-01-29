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

std::shared_ptr<Prover> make_prover(Engine e,
                                    const Property & p,
                                    const TransitionSystem & ts,
                                    const smt::SmtSolver & slv,
                                    PonoOptions opts = PonoOptions());

std::shared_ptr<Prover> make_ceg_proph_prover(Engine e,
                                              const Property & p,
                                              const TransitionSystem & ts,
                                              const smt::SmtSolver & slv,
                                              PonoOptions opts = PonoOptions());

std::shared_ptr<Prover> make_cegar_values_prover(
    Engine e,
    const Property & p,
    const TransitionSystem & ts,
    const smt::SmtSolver & slv,
    PonoOptions opts = PonoOptions());

std::shared_ptr<Prover> make_cegar_bv_arith_prover(
    Engine e,
    const Property & p,
    const TransitionSystem & ts,
    const smt::SmtSolver & slv,
    PonoOptions opts = PonoOptions());

}  // namespace pono
