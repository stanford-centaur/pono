/*********************                                                  */
/*! \file ic3base.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan, Florian Lonsing
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract base class implementation of IC3 parameterized by
**        the unit used in frames, pre-image computation, and inductive
**        and predecessor generalization techniques.
**
**/
#pragma once

#include "engines/prover.h"

namespace pono {

template <class Unit>
class IC3Base : public Prover
{
 public:
  IC3Base(Property & p, smt::SolverEnum se);
  IC3Base(Property & p, const smt::SmtSolver & s);
  IC3Base(const PonoOptions & opt, Property & p, smt::SolverEnum se);
  IC3Base(const PonoOptions & opt, Property & p, const smt::SmtSolver & s);

  typedef Prover super;

  void initialize() override;

 protected:
  ///< the frames data structure.
  ///< a vector of the given Unit template
  ///< which changes depending on the implementation
  std::vector<std::vector<Unit>> frames_;
};

}  // namespace pono
