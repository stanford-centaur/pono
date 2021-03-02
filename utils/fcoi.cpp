/*********************                                                  */
/*! \file coi.cpp
** \verbatim
** Top contributors (to current version):
**   Florian Lonsing, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Class for performing cone of influence reduction
**
**
**/

#include "utils/fcoi.h"

#include "assert.h"
#include "smt-switch/utils.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

FunctionalConeOfInfluence::FunctionalConeOfInfluence(
    const TransitionSystem & ts, int verbosity)
    : ts_(ts), verbosity_(verbosity), local_logger_(verbosity_)
{
  if (!ts_.is_functional()) {
    throw PonoException(
        "Temporary restriction: cone-of-influence analysis "
        "currently supported for functional transition systems only.");
  }
}

/* Main COI function. */
void FunctionalConeOfInfluence::compute_coi(const TermVec & terms)
{
  // clear all the data structures from a previous call
  clear();

  assert(coi_visited_terms_.empty());
  if (verbosity_ >= 3) print_coi_info(terms);

  UnorderedTermSet new_coi_state_vars;
  UnorderedTermSet new_coi_input_vars;

  /* Traverse terms and collect all state/input variables. */
  local_logger_.log(1, "COI analysis");

  for (auto t : terms) {
    compute_term_coi(t, new_coi_state_vars, new_coi_input_vars);
  }

  assert(statevars_in_coi_.empty());
  assert(inputvars_in_coi_.empty());

  /* Add state/input variables found in bad-state term to global collections. */
  for (auto sv : new_coi_state_vars) statevars_in_coi_.insert(sv);
  for (auto sv : new_coi_input_vars) inputvars_in_coi_.insert(sv);

  /* Traverse constraints and collect all state/input variables. */
  local_logger_.log(1, "COI analysis: constraints");
  compute_coi_trans_constraints();

  /* Traverse next-state functions of state-variables that were
     already collected. The loop breaks when no new state/input
     variables were found. Every term is visited at most once during
     the whole COI analysis. */
  unsigned int num_statevars;
  unsigned int num_inputvars;
  unsigned int iterations = 0;
  do {
    iterations++;
    num_statevars = statevars_in_coi_.size();
    num_inputvars = inputvars_in_coi_.size();

    local_logger_.log(
        1, "COI analysis: next-state functions, iteration {}", iterations);
    compute_coi_next_state_funcs();

  } while (statevars_in_coi_.size() != num_statevars
           || inputvars_in_coi_.size() != num_inputvars);

  /* TODO/NOTE: we do NOT traverse 'init_' constraint to search for
     new state variables. The initial constraint term can have any
     arbitrary structure and hence may be difficult to analyze
     precisely. */

  if (verbosity_ >= 3) {
    local_logger_.log(3, "COI analysis completed");
    for (auto var : statevars_in_coi_)
      local_logger_.log(3, "  - found COI statevar {}", var);
    for (auto var : inputvars_in_coi_)
      local_logger_.log(3, "  - found COI inputvar {}", var);

    local_logger_.log(3, "Original system had:");
    for (auto var : ts_.statevars())
      local_logger_.log(3, "  - statevar {}", var);
    for (auto var : ts_.inputvars())
      local_logger_.log(3, "  - inputvar {}", var);
  }
}

void FunctionalConeOfInfluence::clear()
{
  statevars_in_coi_.clear();
  inputvars_in_coi_.clear();
  coi_visited_terms_.clear();
}

/* For debugging only. */
void FunctionalConeOfInfluence::print_term_dfs(const Term & term)
{
  UnorderedTermSet visited_terms;
  TermVec open_terms;
  open_terms.push_back(term);
  Term cur;

  while (!open_terms.empty()) {
    cur = open_terms.back();
    open_terms.pop_back();

    if (visited_terms.find(cur) == visited_terms.end()) {
      // cache 'cur' and push its children
      visited_terms.insert(cur);

      cout << "  visiting term: " << cur << "\n";
      if (cur->is_symbol()) cout << "    ..is symbol\n";

      for (auto child : cur) {
        cout << "    pushing child: " << child << "\n";
        open_terms.push_back(child);
      }
    }
  }
}

/* For debugging only. */
void FunctionalConeOfInfluence::print_coi_info(const TermVec & terms)
{
  cout << "TEST PRINT COI\n";
  cout << "terms to keep:" << endl;
  for (auto t : terms) {
    cout << t << endl;
    print_term_dfs(t);
  }

  cout << "init_ term: " << ts_.init() << "\n";
  print_term_dfs(ts_.init());

  cout << "trans_ term: " << ts_.trans() << "\n";
  print_term_dfs(ts_.trans());

  cout << "input vars: \n";
  for (auto inputvar : ts_.inputvars()) cout << "  " << inputvar << "\n";

  cout << "state vars: \n";
  for (auto statevar : ts_.statevars()) cout << "  " << statevar << "\n";

  cout << "constraints: \n";
  for (auto constr : ts_.constraints()) cout << "  " << constr.first << "\n";
}

/* Add 'term' to 'set' if it does not already appear there. */
/* Traverse 'term', collect state/input variables, and add them to
   global sets 'new_coi_state_vars' and 'new_coi_input_vars'. */
void FunctionalConeOfInfluence::compute_term_coi(
    const Term & term,
    UnorderedTermSet & new_coi_state_vars,
    UnorderedTermSet & new_coi_input_vars)
{
  assert(term != NULL);
  TermVec open_terms;
  open_terms.push_back(term);
  Term cur;

  /* Depth-first search of term structure 'term'. */
  while (!open_terms.empty()) {
    cur = open_terms.back();
    open_terms.pop_back();

    /* Use global set 'coi_visited_terms' here to avoid visiting terms
       multiple times when we call this function on different terms. */
    if (coi_visited_terms_.find(cur) == coi_visited_terms_.end()) {
      /* Cache 'cur' and push its children. */
      coi_visited_terms_.insert(cur);
      local_logger_.log(3, "  visiting COI term: {}", cur);
      if (cur->is_symbol()) {
        local_logger_.log(3, "    ..is symbol");
        if (ts_.statevars().find(cur) != ts_.statevars().end()) {
          assert(ts_.inputvars().find(cur) == ts_.inputvars().end());
          assert(!ts_.is_next_var(cur));
          local_logger_.log(3, "collect COI statevar {}", cur);
          new_coi_state_vars.insert(cur);
        } else if (ts_.inputvars().find(cur) != ts_.inputvars().end()) {
          assert(!ts_.is_curr_var(cur));
          assert(!ts_.is_next_var(cur));
          local_logger_.log(3, "collect COI inputvar {}", cur);
          new_coi_input_vars.insert(cur);
        }
      }

      for (auto child : cur) {
        local_logger_.log(3, "    pushing child: {}", child);
        open_terms.push_back(child);
      }
    }
  }
}

/* Collect state/input variables that appear in next-state functions
   of state-variables that were already collected. */
void FunctionalConeOfInfluence::compute_coi_next_state_funcs()
{
  UnorderedTermSet new_coi_state_vars;
  UnorderedTermSet new_coi_input_vars;
  /* Seed the search using state-variables that were collected already. */
  TermVec unprocessed_state_vars;
  for (auto state_var : statevars_in_coi_)
    unprocessed_state_vars.push_back(state_var);

  while (!unprocessed_state_vars.empty()) {
    Term state_var = unprocessed_state_vars.back();
    unprocessed_state_vars.pop_back();
    assert(ts_.is_curr_var(state_var));
    const smt::UnorderedTermMap & state_updates = ts_.state_updates();
    Term next_func = NULL;
    /* May find state variables without next-function. */
    auto elem = ts_.state_updates().find(state_var);
    if (elem != ts_.state_updates().end()) next_func = elem->second;

    assert(new_coi_state_vars.empty());
    assert(new_coi_input_vars.empty());

    if (next_func != NULL)
      compute_term_coi(next_func, new_coi_state_vars, new_coi_input_vars);

    /* Add newly found state/input variables to global sets
       'statevars_in_coi_' and 'inputvars_in_coi_'. */
    for (auto sv : new_coi_state_vars) {
      if (statevars_in_coi_.find(sv) == statevars_in_coi_.end()) {
        statevars_in_coi_.insert(sv);
        unprocessed_state_vars.push_back(sv);
      }
    }
    for (auto sv : new_coi_input_vars) {
      inputvars_in_coi_.insert(sv);
    }

    new_coi_state_vars.clear();
    new_coi_input_vars.clear();
  }
}

/* Collect state/input variables that appear in constraints that were
   added to the transition system. */
void FunctionalConeOfInfluence::compute_coi_trans_constraints()
{
  UnorderedTermSet new_coi_state_vars;
  UnorderedTermSet new_coi_input_vars;

  for (const auto & e : ts_.constraints()) {
    local_logger_.log(3, "  trans constraints--constr: {}", e.first);
    compute_term_coi(e.first, new_coi_state_vars, new_coi_input_vars);
  }

  /* Add newly collected variables to global collections. */
  for (auto sv : new_coi_state_vars)
    if (statevars_in_coi_.find(sv) == statevars_in_coi_.end())
      statevars_in_coi_.insert(sv);

  for (auto sv : new_coi_input_vars)
    if (inputvars_in_coi_.find(sv) == inputvars_in_coi_.end())
      inputvars_in_coi_.insert(sv);
}

}  // namespace pono
