/*********************                                                  */
/*! \file ic3bits.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Bit-level IC3 implementation that splits bitvector variables
**        into the individual bits for bit-level cubes/clauses
**        However, the transition system itself still uses bitvectors
**/

#pragma once

#include "engines/ic3.h"

namespace pono {

class IC3Bits : public IC3
{
 public:
  // itp_se is the SolverEnum for the interpolator

  IC3Bits(const SafetyProperty & p,
          const TransitionSystem & ts,
          const smt::SmtSolver & s,
          PonoOptions opt = PonoOptions());

  virtual ~IC3Bits() {}

  typedef IC3 super;

  void initialize() override;

 protected:
  smt::TermVec state_bits_;  ///< boolean variables + bit-vector variables
                             ///< split into individual bits

  // virtual method overrides

  IC3Formula get_model_ic3formula() const override;

  bool ic3formula_check_valid(const IC3Formula & u) const override;

  void check_ts() const override;
};

}  // namespace pono
