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
#include <unordered_map>
#include <vector>

#include "smt-switch/utils.h"
#include "utils/logger.h"
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

/** A conjunct next_var = update_expression of a transition relation, read as
 *  a definition of the update of the next-state variable next_var. It
 *  depends on the next-state variables in update_expression, and gives the
 *  update once they all have one. A conjunct x' = y' is read both ways, as
 *  two definitions. */
struct StateUpdateDefinition
{
  size_t conjunct_id;      ///< the index of the conjunct
  Term next_var;           ///< the next-state variable on one side
  Term update_expression;  ///< the other side
  TermVec dependencies;    ///< the next-state variables in update_expression
  size_t num_unresolved;   ///< how many of the dependencies are unresolved,
                           ///< at most dependencies.size()
};

/** Takes the updates that conjuncts of a transition relation define. A
 *  conjunct x' = f(X, I, X'), with the sides in either order, defines the
 *  update x' := f, which must not refer to next-state variables. So each
 *  next-state variable y' in f, on which the definition depends, is replaced
 *  by the update of y', which resolves the dependency. E.g., x' = f(X, I, y')
 *  and y' = g(X, I) resolve to x' := f(X, I, g(X, I)).
 *
 *  A definition that depends on no next-state variable is ready: it gives
 *  its variable the update, unless another definition gave it one before.
 *  Each definition that depends on this variable then has one unresolved
 *  dependency less, and is ready once it has none left. So each conjunct is
 *  visited once, whatever the order of the conjuncts, and each update is
 *  substituted once.
 *  Definitions that depend on each other in a cycle, or on a next-state
 *  variable that no conjunct defines, stay unresolved.
 *
 *  @param ts the transition system
 *  @param conjuncts conjuncts of its transition relation
 *  @param rest gets the conjuncts that define no update
 *  @return the updates, by next-state variable
 */
UnorderedTermMap take_updates(const TransitionSystem & ts,
                              const TermVec & conjuncts,
                              TermVec & rest)
{
  vector<StateUpdateDefinition> definitions;
  // the definitions that depend on each next-state variable
  unordered_map<Term, vector<size_t>> dependents;
  // the definitions that got ready, in that order
  vector<size_t> ready;
  for (size_t i = 0; i < conjuncts.size(); ++i) {
    Term conjunct = conjuncts[i];
    // Currently, only an equality between two terms can define an update. Any
    // other conjunct goes to `rest`, as would a chain a = b = c that a solver
    // kept as one term instead of a = b and b = c.
    TermVec sides;
    if (conjunct->get_op().prim_op == Equal) {
      for (const Term & side : conjunct) {
        sides.push_back(side);
      }
    }
    // Each side that is a next-state variable gives a definition, so x' = y'
    // gives two, x' := y' and y' := x', of which at most one gives an update.
    for (size_t j = 0; sides.size() == 2 && j < 2; ++j) {
      if (!ts.is_next_var(sides[j])) {
        continue;
      }
      UnorderedTermSet symbols;
      get_free_symbolic_consts(sides[1 - j], symbols);
      TermVec dependencies;
      for (const Term & symbol : symbols) {
        if (ts.is_next_var(symbol)) {
          dependencies.push_back(symbol);
        }
      }
      const size_t d = definitions.size();
      for (const Term & next : dependencies) {
        dependents[next].push_back(d);
      }
      if (dependencies.empty()) {
        ready.push_back(d);
      }
      definitions.push_back(
          { i, sides[j], sides[1 - j], dependencies, dependencies.size() });
    }
  }

  const SmtSolver & solver = ts.solver();
  UnorderedTermMap updates;
  // whether each conjunct gave an update, as the others go to `rest`
  vector<bool> gave_update(conjuncts.size(), false);
  for (size_t k = 0; k < ready.size(); ++k) {
    const StateUpdateDefinition & definition = definitions[ready[k]];
    if (updates.find(definition.next_var) != updates.end()) {
      continue;  // a second equality for the same variable goes to `rest`
    }
    // only the updates of the dependencies, as a solver may go through the
    // whole map on each substitution (e.g., Bitwuzla)
    UnorderedTermMap substitution;
    for (const Term & next : definition.dependencies) {
      substitution[next] = updates.at(next);
    }
    updates[definition.next_var] =
        substitution.empty()
            ? definition.update_expression
            : solver->substitute(definition.update_expression, substitution);
    gave_update[definition.conjunct_id] = true;
    // The update resolves the variable for each definition that depends on
    // it: that definition has one unresolved dependency less, and is ready
    // once it has none left.
    const auto it = dependents.find(definition.next_var);
    if (it != dependents.end()) {
      for (size_t d : it->second) {
        if (--definitions[d].num_unresolved == 0) {
          ready.push_back(d);
        }
      }
    }
  }

  for (size_t i = 0; i < conjuncts.size(); ++i) {
    if (!gave_update[i]) {
      rest.push_back(conjuncts[i]);
    }
  }
  return updates;
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

  // Next-state variables occur in the transition relation only. The state
  // updates a relational system records are not used, as setting its
  // transition relation anew does not drop them.
  TermVec rest;  // the conjuncts that define no update
  UnorderedTermMap next_values =
      take_updates(ts, trans_without_constraints(ts), rest);

  // In `rest`, an input chooses each next value that no conjunct defined. A
  // state variable whose next-state variable does not occur in `rest` either
  // takes any value in the next state that the constraints allow, as one
  // without an update does.
  Term rest_term;
  UnorderedTermSet rest_symbols;
  if (!rest.empty()) {
    rest_term = rest.front();
    for (size_t i = 1; i < rest.size(); ++i) {
      rest_term = fts.make_term(And, rest_term, rest[i]);
    }
    get_free_symbolic_consts(rest_term, rest_symbols);
  }
  for (const auto & sv : ts.statevars()) {
    const Term next = ts.next(sv);
    if (next_values.find(next) != next_values.end()) {
      fts.assign_next(sv, next_values.at(next));
    } else if (rest_symbols.find(next) != rest_symbols.end()) {
      // x_in, which chooses the next value of x
      const Term choice =
          allow_new_inputvars
              ? fts.make_generated_inputvar("in", sv, sv->get_sort())
              : fts.make_generated_statevar("in", sv, sv->get_sort());
      fts.assign_next(sv, choice);
      next_values[next] = choice;
    }
  }
  // So far, the only new variables are the x_in.
  logger.log(1,
             "Functionalized: {} of {} state variables have an update, "
             "{} of them chosen by {}",
             fts.state_updates().size(),
             ts.statevars().size(),
             fts.inputvars().size() + fts.statevars().size()
                 - ts.inputvars().size() - ts.statevars().size(),
             allow_new_inputvars ? "an input" : "a state variable");

  if (rest_term) {
    const Term valid =
        fts.make_generated_statevar("valid", fts.make_sort(BOOL));
    const Term step = fts.solver()->substitute(rest_term, next_values);
    fts.assign_next(valid, fts.make_term(And, valid, step));
    // Also starts valid off true.
    fts.add_constraint(valid);
  }

  // the constraints that `rest` leaves out
  for (const auto & [constraint, to_init_and_next] : ts.constraints()) {
    fts.add_constraint(constraint, to_init_and_next);
  }

  return fts;
}

}  // namespace pono
