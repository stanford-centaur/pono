/*********************                                                        */
/*! \file nondet_fault_injector.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
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

#include "modifiers/nondet_fault_injector.h"

#include "utils/logger.h"

using namespace std;
using namespace smt;

namespace pono {

void NonDetFaultInjector::create_fault_vals()
{
  Term faultval;
  Term st;
  for (auto elem : fts_.state_updates()) {
    st = elem.first;
    faultval = faulty_fts_.make_inputvar("faultval_" + st->to_string(),
                                         st->get_sort());
    state2faultval_[st] = faultval;
    faultval2state_[faultval] = st;
  }
}

}  // namespace pono
