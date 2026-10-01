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
#include "smt-switch/smt.h"

namespace pono {

namespace moxi {
struct Check;
struct Query;
struct SortInfo;
struct System;
class Reader;
}  // namespace moxi

/** Encodes the query of a check-system command of a MoXI file as a safety
 *  property of a relational transition system.
 *
 *  A query brings its own assumptions and initiality condition into the
 *  transition system, which therefore holds a single query: the command must
 *  have exactly one.
 *
 *  The query's system is flattened into the transition system: each
 *  subsystem instance contributes its own copy of its local variables, named
 *  after the instance, e.g. inst.var. The variables take the names that the
 *  check-system command gives them.
 *
 *  The property holds exactly if the query is unsatisfiable, i.e., if no
 *  finite trace starting in an initial state (or in a state satisfying the
 *  query's :current formula, which replaces the initial condition) and
 *  keeping the system's invariant and the query's assumptions reaches each
 *  of its reachability conditions, possibly at different steps.
 *
 *  Unlike the n-satisfiability that the MoXI description defines, the trace
 *  does not need to extend by another step past its last state; this is the
 *  usual semantics of reachability, which e.g. MoXIchecker implements too.
 *  The two only differ if the transition relation deadlocks.
 *
 *  Queries with fairness conditions ask for infinite traces and are not
 *  supported yet.
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
   */
  MoxiEncoder(const std::string & filename,
              RelationalTransitionSystem & rts,
              std::size_t check_idx = 0);

  ~MoxiEncoder();

  /** @return the property, which fails exactly if the query is satisfiable.
   *  It has next-state variables if a reachability condition does, which
   *  pono accounts for by monitoring it. */
  const smt::Term & prop() const { return prop_; }

  /** @return the name of the encoded query */
  const std::string & query_name() const { return query_name_; }

 private:
  /** A variable of the flattened system. */
  struct FlatVar
  {
    std::string name;
    const moxi::SortInfo * sort;
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
                          const moxi::SortInfo & sort,
                          bool is_state);

  RelationalTransitionSystem & rts_;
  std::unique_ptr<moxi::Reader> reader_;
  std::string query_name_;
  smt::Term prop_;

  // the flattened system
  std::deque<FlatVar> vars_;  ///< a deque keeps pointers valid
  std::vector<Instance> instances_;
  std::unordered_set<std::string> taken_names_;
};

}  // namespace pono
