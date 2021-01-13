/*********************                                                  */
/*! \file ic3.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Bit-level IC3 implementation using the IC3Base abstract base class
**/

#pragma once

#include "engines/ic3base.h"

namespace pono {

class IC3 : public IC3Base
{
 public:
  IC3(const Property & p, const TransitionSystem & ts,
      const smt::SmtSolver & s, PonoOptions opt = PonoOptions());

  virtual ~IC3() {}

  typedef IC3Base super;

 protected:
  // pure virtual method implementations

  IC3Formula get_model_ic3formula() const override;

  bool ic3formula_check_valid(const IC3Formula & u) const override;

  void predecessor_generalization(size_t i,
                                  const smt::Term & c,
                                  IC3Formula & pred) override;

  void check_ts() const override;

};

}  // namespace pono
