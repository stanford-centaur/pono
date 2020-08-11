/*********************                                                  */
/*! \file prophecy_refiner.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Class for adding prophecy variables that predict a value some
**        constant number of steps before a property violation.
**/

#include "modifiers/prophecy_refiner.h"

using namespace smt;
using namespace std;

namespace pono {

ProphecyRefiner::ProphecyRefiner(TransitionSystem & ts)
    : ts_(ts), solver_(ts_.solver()), hr_(ts_)
{
}

pair<Term, Term> ProphecyRefiner::get_proph(const Term & target,
                                            size_t delay,
                                            const Term & prop)
{
  // first use history variables to delay target
  Term hist_var = hr_.get_hist(target, delay);

  // now add a prophecy variable which targets that history variable
  string name = "proph_" + target->to_string() + "_" + std::to_string(delay);
  Term proph_var = ts_.make_statevar(name, target->get_sort());
  // make it frozen
  ts_.assign_next(proph_var, proph_var);

  // prophecy variable updates property using an implication
  // (proph = hist) -> prop
  Term antecedent = solver_->make_term(Equal, proph_var, hist_var);
  Term updated_prop = solver_->make_term(Implies, antecedent, prop);
  return { proph_var, updated_prop };
}

}  // namespace pono
