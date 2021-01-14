/*********************                                                        */
/*! \file cegar_values.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief A simple CEGAR loop that abstracts values with frozen variables
**        and refines by constraining the variable to the value again
**
**/

#include "engines/cegar_values.h"

#include "utils/exceptions.h"

using namespace smt;
using namespace std;

namespace pono {

// TODO add a value abstractor
//      make sure not to introduce nonlinearities
//      implement generic backend
//      then specialize for IC3IA

template <class Prover_T>
CegarValues<Prover_T>::CegarValues(const Property & p,
                                   const TransitionSystem & ts,
                                   const smt::SmtSolver & solver,
                                   PonoOptions opt)
    : super(p, ts, solver, opt)
{
}

template <class Prover_T>
ProverResult CegarValues<Prover_T>::check_until(int k)
{
  throw PonoException("NYI");
}

template <class Prover_T>
void CegarValues<Prover_T>::initialize()
{
  throw PonoException("NYI");
}

template <class Prover_T>
void CegarValues<Prover_T>::cegar_abstract()
{
  throw PonoException("NYI");
}

template <class Prover_T>
bool CegarValues<Prover_T>::cegar_refine()
{
  throw PonoException("NYI");
}

}  // namespace pono
