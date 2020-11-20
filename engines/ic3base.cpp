/*********************                                                  */
/*! \file ic3base.cpp
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

#include <type_traits>

#include "engines/ic3base.h"

using namespace smt;
using namespace std;

namespace pono {

IC3Base::IC3Base(Property & p, smt::SolverEnum se, IC3UnitCreator ic)
    : super(p, se), mk_unit(ic)
{
  initialize();
}

IC3Base::IC3Base(Property & p, const smt::SmtSolver & s, IC3UnitCreator ic)
    : super(p, s), mk_unit(ic)
{
  initialize();
}

IC3Base::IC3Base(const PonoOptions & opt,
                 Property & p,
                 smt::SolverEnum se,
                 IC3UnitCreator ic)
    : super(opt, p, se), mk_unit(ic)
{
  initialize();
}
IC3Base::IC3Base(const PonoOptions & opt,
                 Property & p,
                 const smt::SmtSolver & s,
                 IC3UnitCreator ic)
    : super(opt, p, s), mk_unit(ic)
{
  initialize();
}

void IC3Base::initialize()
{
  // TODO: fix initialize. Not sure it makes sense to call it here instead of in
  // prover constructor
  //       need to think about multiple levels of inheritance and where the
  //       responsibility for calling initialize belongs
  super::initialize();
}

}  // namespace pono
