/*********************                                                        */
/*! \file fault_injector.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the cosa2 project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Modifies a transition system such that faults can occur on state
**        variables at each timestep. This assumes a functional model, and
**        inserts a mux and a boolean FAULT variable for each state element.
**        If the fault variable is true, it can get a non-deterministic value
**
**/

#include "fault_injector.h"
#include "logger.h"

using namespace std;
using namespace smt;

namespace cosa
{

void FaultInjector::do_fault_injection()
{
  const UnorderedTermMap & state_updates = fts_.state_updates();
  SmtSolver & solver = faulty_fts_.solver();

  Sort boolsort = solver->make_sort(BOOL);
  Term s;
  Term faultsig;
  Term faultsel;
  Term faultval;
  for (auto elem : state_updates)
  {
    s = elem.first;

    faultsel = faulty_fts_.make_input("faultsel_" + s->to_string(), boolsort);
    state2faultsel_[s] = faultsel;
    faultsel2state_[faultsel] = s;

    faultsig = faulty_fts_.make_state("FAULT_" + s->to_string(), boolsort);
    state2faultsig_[s] = faultsig;
    faultsig2state_[faultsig] = s;

    faulty_fts_.constrain_init(solver->make_term(Equal, faultsig, solver->make_term(false)));
    faulty_fts_.set_next(faultsig, solver->make_term(Or, faultsig, faultsel));

    faultval = faulty_fts_.make_input("faultval_" + s->to_string(), s->get_sort());
    state2faultval_[s] = faultval;
    faultval2state_[faultval] = s;

    faulty_fts_.set_next(s, solver->make_term(Ite,
                                              faultsig,
                                              faultval,
                                              elem.second));
    fault_sigs_.push_back(faultsig);
  }
}

}
