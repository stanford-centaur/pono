/*********************                                                        */
/*! \file ceg_prophecy.cpp
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

#include "engines/ceg_prophecy.h"

using namespace smt;
using namespace std;

namespace pono {

CegProphecy::CegProphecy(const Property & p, smt::SolverEnum se)
    : super(p, se),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, unroller_)
{
  super::initialize();
}

CegProphecy::CegProphecy(const Property & p, const SmtSolver & solver)
    : super(p, solver),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, unroller_)
{
  super::initialize();
}

CegProphecy::CegProphecy(const PonoOptions & opt,
                         const Property & p,
                         smt::SolverEnum se)
    : super(opt, p, se),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, unroller_)
{
  super::initialize();
}

CegProphecy::CegProphecy(const PonoOptions & opt,
                         const Property & p,
                         const smt::SmtSolver & solver)
    : super(opt, p, solver),
      conc_ts_(p.transition_system()),
      solver_(conc_ts_.solver()),
      abs_ts_(solver_),
      aa_(conc_ts_, abs_ts_, true),
      aae_(p, aa_, unroller_)
{
  super::initialize();
}

ProverResult CegProphecy::check_until(int k) { throw PonoException("NYI"); }

void CegProphecy::abstract() { throw PonoException("NYI"); }

void CegProphecy::refine() { throw PonoException("NYI"); }

}  // namespace pono
