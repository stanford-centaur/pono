/*********************                                                  */
/*! \file functionalize.cpp
** \verbatim
** Top contributors (to current version):
**   Po-Chun Chien
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Turns a relational transition system into a functional one
**
**
**/

#include "modifiers/functionalize.h"

#include <algorithm>

#include "smt-switch/utils.h"
#include "utils/ts_manipulation.h"

using namespace smt;
using namespace std;

namespace pono {

namespace {

/** Returns the conjuncts of the transition relation of a system, except for
 *  true and the conjuncts that its constraints added. The functional system
 *  keeps the constraints as its own, because engines see them apart from the
 *  transition relation, e.g. IC3SA in its abstraction.
 *  @param ts the transition system
 *  @return the conjuncts
 */
TermVec trans_without_constraints(const TransitionSystem & ts)
{
  // as TransitionSystem::add_constraint adds them: a constraint, which has
  // no next-state variables, and its next-state copy if it is to hold in the
  // next state and refers to state variables only
  TermVec added{ ts.solver()->make_term(true) };
  for (const auto & [constraint, to_init_and_next] : ts.constraints()) {
    conjunctive_partition(constraint, added, true);
    if (to_init_and_next && ts.only_curr(constraint)) {
      conjunctive_partition(ts.next(constraint), added, true);
    }
  }
  const UnorderedTermSet skip(added.begin(), added.end());
  TermVec conjuncts;
  conjunctive_partition(ts.trans(), conjuncts, true);
  conjuncts.erase(remove_if(conjuncts.begin(),
                            conjuncts.end(),
                            [&skip](const Term & conjunct) {
                              return skip.find(conjunct) != skip.end();
                            }),
                  conjuncts.end());
  return conjuncts;
}

}  // namespace

TransitionSystem functionalize(const TransitionSystem & ts,
                               bool allow_new_inputvars)
{
  if (ts.is_functional()) {
    return ts;
  }

  TransitionSystem fts = create_fresh_ts(true, ts.solver());
  for (const auto & sv : ts.statevars()) {
    fts.add_statevar(sv, ts.next(sv));
  }
  for (const auto & iv : ts.inputvars()) {
    fts.add_inputvar(iv);
  }
  for (const auto & [name, term] : ts.named_terms()) {
    fts.name_term(name, term);
  }
  fts.set_init(ts.init());

  // Next-state variables occur in the transition relation only. A state
  // variable whose next-state variable occurs in T only in a constraint takes
  // any value in the next state that the constraints allow, as one without
  // an update does. The state updates a relational system records are not
  // used, as setting its transition relation anew does not drop them.
  Term relation = fts.make_term(true);
  for (const Term & conjunct : trans_without_constraints(ts)) {
    relation = fts.make_term(And, relation, conjunct);
  }
  UnorderedTermSet relation_symbols;
  get_free_symbolic_consts(relation, relation_symbols);
  UnorderedTermMap next_values;
  for (const auto & sv : ts.statevars()) {
    const Term next = ts.next(sv);
    if (relation_symbols.find(next) == relation_symbols.end()) {
      continue;
    }
    // x_in, which chooses the next value of x
    const Term choice =
        allow_new_inputvars
            ? fts.make_generated_inputvar("in", sv, sv->get_sort())
            : fts.make_generated_statevar("in", sv, sv->get_sort());
    fts.assign_next(sv, choice);
    next_values[next] = choice;
  }

  const Term valid = fts.make_generated_statevar("valid", fts.make_sort(BOOL));
  const Term step = fts.solver()->substitute(relation, next_values);
  fts.assign_next(valid, fts.make_term(And, valid, step));
  // Also starts valid off true.
  fts.add_constraint(valid);

  // the constraints that the step leaves out
  for (const auto & [constraint, to_init_and_next] : ts.constraints()) {
    fts.add_constraint(constraint, to_init_and_next);
  }

  return fts;
}

}  // namespace pono
