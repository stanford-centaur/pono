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

  faulty_fts_.set_trans(solver->make_term(true));

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
    faultsel2val_[faultsel] = s;

    faultsig = faulty_fts_.make_state("FAULT_" + s->to_string(), boolsort);
    state2faultsig_[s] = faultsig;
    faultsig2val_[faultsig] = s;
    faultsig2name_[faultsig] = s->to_string();

    faulty_fts_.constrain_init(solver->make_term(Equal, faultsig, solver->make_term(false)));
    faulty_fts_.set_next(faultsig, solver->make_term(Or, faultsig, faultsel));

    faultval = faulty_fts_.make_input("faultval_" + s->to_string(), s->get_sort());
    state2faultval_[s] = faultval;
    faultval2val_[faultval] = s;

    faulty_fts_.set_next(s, solver->make_term(Ite,
                                              faultsig,
                                              faultval,
                                              elem.second));
    fault_sigs_.push_back(faultsig);
  }

  Term val;
  string name;
  for (auto elem : fts_.named_terms())
  {
    name = elem.first;
    val = elem.second;

    faultsel = faulty_fts_.make_state("faultsel_" + name, boolsort);

    faultsig = faulty_fts_.make_state("FAULT_" + name, boolsort);
    faultsig2name_[faultsig] = name;
    faultsig2val_[faultsig] = val;

    faulty_fts_.constrain_init(solver->make_term(Equal, faultsig, solver->make_term(false)));

    Term update = solver->make_term(Or, faulty_fts_.next(faultsel),
                                    solver->make_term(Or, faultsig, faultsel));
    faulty_fts_.constrain_trans(solver->make_term(Equal,
                                                  faulty_fts_.next(faultsig),
                                                  update));
    // faulty_fts_.set_next(faultsig, solver->make_term(Or, faultsig, faultsel));

    faultval = faulty_fts_.make_input("faultval_" + name, val->get_sort());

    faulty_fts_.remove_name(name);
    faulty_fts_.name_term(name,
                          solver->make_term(Ite,
                                            faultsel,
                                            faultval,
                                            val));
    fault_sigs_.push_back(faultsig);
  }
}

}
