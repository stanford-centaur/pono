/*********************                                                        */
/*! \file moxi_encoder.h
** \verbatim
** Top contributors (to current version):
**   Po-Chun Chien
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** (in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Frontend for the Model Exchange Interlingua (MoXI).
**        See https://doi.org/10.1007/978-3-031-65627-9_10 and
**        https://github.com/ModelChecker/IL for more information.
**
**/

#pragma once

#include <cstddef>
#include <deque>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "core/rts.h"
#include "options/options.h"
#include "smt-switch/smt.h"

namespace pono {

namespace moxi {
struct Check;
struct MoxiSort;
struct Query;
struct System;
class Reader;
}  // namespace moxi

/** Encodes the query of a check-system command of a MoXI file as a safety
 *  property of a relational transition system.
 *
 *  A query brings its own assumptions and initiality condition into the
 *  transition system, which therefore holds a single query: the command must
 *  have exactly one, from a :query attribute. The :queries attribute is not
 *  supported, as its queries must share the values of the declared constants
 *  and functions, which checking them one at a time cannot ensure.
 *
 *  The query's system is flattened into the transition system: each
 *  subsystem instance contributes its own copy of its local variables, named
 *  after the instance, e.g. inst.var. The variables take the names that the
 *  check-system command gives them.
 *
 *  The property holds exactly if the query is unsatisfiable, i.e., not
 *  n-satisfiable for any n, as the MoXI description defines it: no trace
 *  starting in an initial state (or in a state satisfying the query's
 *  :current formula, which replaces the initial condition) takes n + 1
 *  steps that keep the system's transition relation and invariant and the
 *  query's assumptions, and meets each of its reachability conditions at one
 *  of these steps, possibly at different ones.
 *
 *  This differs from the usual semantics of reachability at the last state
 *  of the trace, the one that the last step leads to:
 *  - It must exist, so even after the conditions hold, the trace needs a
 *    successor by the transition relation and the assumptions over
 *    next-state inputs. The description requires systems to be
 *    deadlock-free, but the assumptions may still allow no next input, and
 *    a state that only the query's :current formula makes reachable may have
 *    no successor.
 *  - It need not satisfy the invariant or the one-state assumptions,
 *    although a reachability condition over next-state variables refers to
 *    it.
 *  The encoding therefore imposes the invariant and the assumptions only on
 *  the state that a step leaves, and a flag remembers for each reachability
 *  condition that it held at a step. The property fails once all flags are
 *  set, so a counterexample ends in that last state, which may violate the
 *  invariant.
 *
 *  A stutter step goes from a state to the same state. It may be added to a
 *  system to make it deadlock-free, as MoXI requires: the transition
 *  relation becomes a disjunction of the original transition relation and a
 *  stutter term, which is a conjunction of equalities (= v' v), e.g.
 *  (and (= x' x) (= y' y)) for a system with the variables x and y. We call
 *  such an equality a stutter equality, and say that it keeps the variable
 *  v. The stutter term has one for each variable, so that its steps are
 *  stutter steps, which the encoding allows from each state that satisfies
 *  the invariant and the assumptions.
 *
 *  With the experimental option --moxi-drop-stutter, the encoder reverses
 *  this: if the transition relation is a disjunction with a stutter term,
 *  it drops that term. As a stutter step reaches no new state, the property
 *  then fails in a state that satisfies the invariant and the assumptions
 *  and in which each reachability condition holds or held before: from such
 *  a state, the stutter step would complete the trace that the query asks
 *  for. A counterexample therefore ends in this state; the stutter step
 *  that would follow only repeats it.
 *
 *  The stutter term is dropped only for a query whose formulas have no
 *  next-state variables and for a system without subsystems, which would
 *  constrain the stutter step. Of the disjuncts that consist of stutter
 *  equalities, the encoder takes the one that keeps the most variables as
 *  the stutter term if it keeps each variable that becomes a state variable.
 *  It need not keep the variables that become inputs, which take values of
 *  their own in each step anyway. Otherwise the transition relation stays
 *  as it is.
 *
 *  Queries with fairness conditions ask for infinite traces and are not
 *  supported yet.
 *
 *  An enumeration sort is a datatype, or another sort, e.g. Int, whose values
 *  0, 1, ... stand for those of the enumeration. Then a variable of the
 *  enumeration takes only these values, in every state, including the last
 *  one, and arrays indexed by the enumeration are compared at these values
 *  only. Where the solver would still tell apart arrays that are equal in
 *  MoXI, the reader refuses the file.
 */
class MoxiEncoder
{
 public:
  /** Parses a MoXI file and encodes the query of one of its check-system
   *  commands.
   *  @param filename the file to read
   *  @param rts the transition system to encode the query into
   *  @param check_idx the check-system command, counting them in the order
   *         of the file
   *  @param enum_encoding the sort to encode enumeration sorts as
   *  @param drop_stutter whether to drop a stutter term from the transition
   *         relation where the query allows it
   */
  MoxiEncoder(const std::string & filename,
              RelationalTransitionSystem & rts,
              std::size_t check_idx = 0,
              MoxiEnumEncoding enum_encoding = MoxiEnumEncoding::DT,
              bool drop_stutter = false);

  ~MoxiEncoder();

  /** @return the property, which fails exactly if the query is satisfiable.
   *  It refers only to current-state variables: the flags of the
   *  reachability conditions and, with the stutter term dropped, the
   *  conditions, the invariant and the assumptions. */
  const smt::Term & prop() const { return prop_; }

  /** @return the name of the encoded query */
  const std::string & query_name() const { return query_name_; }

  /** @return whether the stutter term was dropped, in which case a
   *  counterexample lacks the stutter step out of its last state */
  bool stutter_dropped() const { return stutter_dropped_; }

 private:
  /** A variable of the flattened system. */
  struct FlatVar
  {
    std::string name;
    const moxi::MoxiSort * sort;  ///< held by the reader
    /** the placeholders of its current and next values in the formulas of
     *  the check-system, and of a constant in all formulas; null for the
     *  local variables of subsystems */
    smt::Term curr_placeholder;
    smt::Term next_placeholder;
    bool frozen;  ///< whether it is a declared constant
    // the references to it, which decide whether it must be a state
    bool in_state = false;      ///< by a formula that holds in each state
    bool in_step = false;       ///< by a formula over transitions
    bool next_in_step = false;  ///< to its next value, by one of the latter
    // its values in the transition system, once it is made
    smt::Term curr = nullptr;
    smt::Term next = nullptr;  ///< remains null for an input or a constant
  };

  /** An instance of a system, i.e. of the check-system's or a subsystem. */
  struct Instance
  {
    const moxi::System * system;
    /** the variables that its variables stand for, in the order
     *  System::variable counts them */
    std::vector<FlatVar *> actuals;
  };

  void encode(const moxi::Check & check, const moxi::Query & query);

  /** Instantiates a system and, recursively, its subsystems, with their own
   *  local variables.
   *  @param system the system to instantiate
   *  @param actuals the variables that its variables stand for, in the order
   *         System::variable counts them
   *  @param prefix the name of the instance, including a trailing dot
   */
  void flatten(const moxi::System & system,
               const std::vector<FlatVar *> & actuals,
               const std::string & prefix);

  /** Creates a variable of the transition system, under the given name
   *  unless a symbol of the solver has it already. */
  smt::Term make_variable(const std::string & name,
                          const smt::Sort & sort,
                          bool is_state);

  RelationalTransitionSystem & rts_;
  std::unique_ptr<moxi::Reader> reader_;
  bool drop_stutter_;
  bool stutter_dropped_ = false;
  std::string query_name_;
  smt::Term prop_;

  // the flattened system
  std::deque<FlatVar> vars_;  ///< a deque keeps pointers valid
  std::vector<Instance> instances_;
  std::unordered_set<std::string> taken_names_;
};

}  // namespace pono
