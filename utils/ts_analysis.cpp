/*********************                                                        */
/*! \file ts_analysis.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Useful functions for analyzing transition systems.
**
**
**/

#include "smt-switch/term_translator.h"

#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/ts_analysis.h"

using namespace smt;

namespace pono {

bool check_invar(const TransitionSystem & ts,
                 const Term & other_prop,
                 const Term & other_invar)
{
  if (!ts.only_curr(other_invar)) {
    logger.log(1, "INVARCHECK: Fail, contains non-current state vars");
    return false;
  }

  // use a fresh solver
  // to avoid issues with a corrupted solver state
  SmtSolver solver = create_solver(ts.solver()->get_solver_enum());
  solver->set_opt("incremental", "true");
  TermTranslator tt(solver);
  Term init = tt.transfer_term(ts.init(), BOOL);
  Term trans = tt.transfer_term(ts.trans(), BOOL);
  Term prop = tt.transfer_term(other_prop, BOOL);
  Term invar = tt.transfer_term(other_invar, BOOL);
  Term next_invar = tt.transfer_term(ts.next(other_invar), BOOL);

  bool pass = true;

  solver->push();
  solver->assert_formula(init);
  solver->assert_formula(solver->make_term(Not, invar));
  Result r = solver->check_sat();
  solver->pop();
  pass &= r.is_unsat();
  logger.log(1, "INVARCHECK: init |= inv...{}", r.is_unsat() ? "OK" : "FAIL");

  solver->push();
  solver->assert_formula(invar);
  solver->assert_formula(trans);
  solver->assert_formula(solver->make_term(Not, next_invar));
  r = solver->check_sat();
  solver->pop();
  pass &= r.is_unsat();
  logger.log(
      1, "INVARCHECK: inv & trans |= inv'...{}", r.is_unsat() ? "OK" : "FAIL");

  solver->push();
  solver->assert_formula(invar);
  solver->assert_formula(solver->make_term(Not, prop));
  r = solver->check_sat();
  solver->pop();
  pass &= r.is_unsat();
  logger.log(1, "INVARCHECK: inv |= prop...{}", r.is_unsat() ? "OK" : "FAIL");

  return pass;
}

}  // namespace pono
