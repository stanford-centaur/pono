/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann, Florian Lonsing
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#include "prover.h"
#include "available_solvers.h"
#include "utils/logger.h"

#include <climits>
#include <cassert>
#include <functional>

using namespace smt;
using namespace std;

namespace pono {

Prover::Prover(Property & p, smt::SolverEnum se)
    : Prover(p, create_solver(se))
{
  solver_->set_opt("incremental", "true");
  solver_->set_opt("produce-models", "true");
}

Prover::Prover(Property & p, const smt::SmtSolver & s)
    : solver_(s),
      to_prover_solver_(s),
      property_(p, to_prover_solver_),
      ts_(property_.transition_system()),
      orig_ts_(p.transition_system()),
      unroller_(ts_, solver_)
{
}

Prover::Prover(const PonoOptions & opt, Property & p, smt::SolverEnum se)
    : Prover(opt, p, create_solver(se))
{
  solver_->set_opt("incremental", "true");
  solver_->set_opt("produce-models", "true");
}

Prover::Prover(const PonoOptions & opt,
               Property & p,
               const smt::SmtSolver & s)
    : solver_(s),
      to_prover_solver_(solver_),
      property_(p, to_prover_solver_),
      ts_(property_.transition_system()),
      orig_ts_(p.transition_system()),
      unroller_(ts_, solver_),
      options_(opt)
{
}

Prover::~Prover() {}

void Prover::initialize()
{
  reached_k_ = -1;
  bad_ = solver_->make_term(smt::PrimOp::Not, property_.prop());
  if (options_.static_coi_) {
    if (!ts_.is_functional())
      throw PonoException("Temporary restriction: cone-of-influence analysis currently "\
                          "supported for functional transition systems only.");
    /* Compute the set of state/input variables related to the
       bad-state property. Based on that information, rebuild the
       transition relation of the transition system. */
    compute_coi();
    orig_num_statevars_ = ts_.statevars().size();
    orig_num_inputvars_ = ts_.inputvars().size();
    ts_.rebuild_trans_based_on_coi(statevars_in_coi_, inputvars_in_coi_);
    assert(statevars_in_coi_.size() == ts_.statevars().size());
    assert(inputvars_in_coi_.size() == ts_.inputvars().size());
    logger.log(
        1,
        "COI analysis completed: {} remaining input variables, {} original",
        inputvars_in_coi_.size(),
        orig_num_inputvars_);
    logger.log(
        1,
        "COI analysis completed: {} remaining state variables, {} original",
        statevars_in_coi_.size(),
        orig_num_statevars_);
  }
}

ProverResult Prover::prove() { return check_until(INT_MAX); }

/* For debugging only. */
void Prover::print_term_dfs(const Term & term)
{
  UnorderedTermSet visited_terms;
  TermVec open_terms;
  open_terms.push_back (term);
  Term cur;

  while (!open_terms.empty()) {
    cur = open_terms.back();
    open_terms.pop_back ();

    if (visited_terms.find(cur) == visited_terms.end()) {
      // cache 'cur' and push its children
      visited_terms.insert(cur);

      cout << "  visiting term: " << cur << "\n";
      if (cur->is_symbol())
        cout << "    ..is symbol\n";

      for (auto child : cur) {
        cout << "    pushing child: " << child << "\n";
        open_terms.push_back(child);
      }
    }
  }
}

/* For debugging only. */
void Prover::print_coi_info()
{
  cout << "TEST PRINT COI\n";
  cout << "bad_ term: " << bad_ << "\n";
  print_term_dfs(bad_);

  cout << "init_ term: " << ts_.init() << "\n";
  print_term_dfs(ts_.init());

  cout << "trans_ term: " << ts_.trans() << "\n";
  print_term_dfs(ts_.trans());

  cout << "input vars: \n";
  for (auto inputvar : ts_.inputvars())
    cout << "  " << inputvar << "\n";

  cout << "state vars: \n";
  for (auto statevar : ts_.statevars())
    cout << "  " << statevar << "\n";

  cout << "constraints: \n";
  for (auto constr : ts_.constraints())
    cout << "  " << constr << "\n";
}

/* Add 'term' to 'set' if it does not already appear there. */
void Prover::collect_coi_term(UnorderedTermSet & set, const Term & term)
{
  if (set.find(term) == set.end())
    set.insert(term);
}

/* Traverse 'term', collect state/input variables, and add them to
   global sets 'new_coi_state_vars' and 'new_coi_input_vars'. */
void Prover::compute_term_coi(const Term & term,
                              UnorderedTermSet & new_coi_state_vars,
                              UnorderedTermSet & new_coi_input_vars)
{
  assert(term != NULL);
  TermVec open_terms;
  open_terms.push_back (term);
  Term cur;

  /* Depth-first search of term structure 'term'. */
  while (!open_terms.empty()) {
    cur = open_terms.back();
    open_terms.pop_back ();

    /* Use global set 'coi_visited_terms' here to avoid visiting terms
       multiple times when we call this function on different terms. */
    if (coi_visited_terms_.find(cur) == coi_visited_terms_.end()) {
      /* Cache 'cur' and push its children. */
      coi_visited_terms_.insert(cur);
      logger.log(3, "  visiting COI term: {}", cur);
      if (cur->is_symbol()) {
        logger.log(3, "    ..is symbol");
        if (ts_.statevars().find(cur) != ts_.statevars().end()) {
          assert(ts_.inputvars().find(cur) == ts_.inputvars().end());
          assert(!ts_.is_next_var(cur));
          logger.log(3, "collect COI statevar {}", cur);
          collect_coi_term(new_coi_state_vars, cur);
        }
        else if (ts_.inputvars().find(cur) != ts_.inputvars().end()) {
          assert(!ts_.is_curr_var(cur));
          assert(!ts_.is_next_var(cur));
          logger.log(3, "collect COI inputvar {}", cur);
          collect_coi_term(new_coi_input_vars, cur);
        }
      }

      for (auto child : cur) {
        logger.log(3, "    pushing child: {}", child);
        open_terms.push_back(child);
      }
    }
  }
}

/* Collect state/input variables that appear in next-state functions
   of state-variables that were already collected. */
void Prover::compute_coi_next_state_funcs()
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
    if (elem != ts_.state_updates().end())
      next_func = elem->second;

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
      collect_coi_term(inputvars_in_coi_, sv);
    }

    new_coi_state_vars.clear();
    new_coi_input_vars.clear();
  }
}

/* Collect state/input variables that appear in constraints that were
   added to the transition system. */
void Prover::compute_coi_trans_constraints()
{
  UnorderedTermSet new_coi_state_vars;
  UnorderedTermSet new_coi_input_vars;

  for (auto constr : ts_.constraints()) {
    logger.log(3, "  trans constraints--constr: {}", constr);
    compute_term_coi(constr, new_coi_state_vars, new_coi_input_vars);
  }

  /* Add newly collected variables to global collections. */
  for (auto sv : new_coi_state_vars)
    if (statevars_in_coi_.find(sv) == statevars_in_coi_.end())
      statevars_in_coi_.insert(sv);

  for (auto sv : new_coi_input_vars)
    if (inputvars_in_coi_.find(sv) == inputvars_in_coi_.end())
      inputvars_in_coi_.insert(sv);
}

/* Main COI function. */
void Prover::compute_coi()
{
  assert (coi_visited_terms_.empty());
  if (options_.verbosity_ >= 3)
    print_coi_info();

  logger.log(1, "Starting static cone-of-influence (COI) analysis:");
  logger.log(1, "  - input variables: {}", ts_.inputvars().size());
  logger.log(1, "  - state variables: {}", ts_.statevars().size());
  logger.log(1, "  - constraints: {}", ts_.constraints().size());

  UnorderedTermSet new_coi_state_vars;
  UnorderedTermSet new_coi_input_vars;

  /* Traverse bad-state property term and collect all state/input variables. */
  logger.log(1, "COI analysis: bad-term");
  compute_term_coi(bad_, new_coi_state_vars, new_coi_input_vars);

  assert (statevars_in_coi_.empty());
  assert (inputvars_in_coi_.empty());

  /* Add state/input variables found in bad-state term to global collections. */
  for (auto sv : new_coi_state_vars)
    statevars_in_coi_.insert(sv);
  for (auto sv : new_coi_input_vars)
    inputvars_in_coi_.insert(sv);

  /* Traverse constraints and collect all state/input variables. */
  logger.log(1, "COI analysis: constraints");
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

    logger.log(1, "COI analysis: next-state functions, iteration {}", iterations);
    compute_coi_next_state_funcs();

  } while (statevars_in_coi_.size() != num_statevars ||
           inputvars_in_coi_.size() != num_inputvars);

  /* TODO/NOTE: we do NOT traverse 'init_' constraint to search for
     new state variables. The initial constraint term can have any
     arbitrary structure and hence may be difficult to analyze
     precisely. */

  if (options_.verbosity_ >= 3) {
    logger.log(3, "COI analysis completed");
    for (auto var : statevars_in_coi_)
      logger.log(3, "  - found COI statevar {}", var);
    for (auto var : inputvars_in_coi_)
      logger.log(3, "  - found COI inputvar {}", var);

    logger.log(3, "Original system had:");
    for (auto var : ts_.statevars())
      logger.log(3, "  - statevar {}", var);
    for (auto var : ts_.inputvars())
      logger.log(3, "  - inputvar {}", var);
  }
}

bool Prover::witness(std::vector<UnorderedTermMap> & out)
{
  if (!witness_.size()) {
    throw PonoException(
        "Recovering witness failed. Make sure that there was "
        "a counterexample and that the engine supports witness generation.");
  }

  function<Term(const Term &, SortKind)> transfer_to_prover_as;
  function<Term(const Term &, SortKind)> transfer_to_orig_ts_as;
  TermTranslator to_orig_ts_solver(orig_ts_.solver());
  if (solver_ == orig_ts_.solver()) {
    // don't need to transfer terms if the solvers are the same
    transfer_to_prover_as = [](const Term & t, SortKind sk) { return t; };
    transfer_to_orig_ts_as = [](const Term & t, SortKind sk) { return t; };
  } else {
    /* TODO: double-check that transferring terms still works as
       intended in this branch when COI is used. */
    if (options_.static_coi_)
      throw PonoException(
          "Temporary restriction: cone-of-influence analysis "
          "currently incompatible with witness generation.");
    // need to add symbols to cache
    UnorderedTermMap & cache = to_orig_ts_solver.get_cache();
    for (auto v : orig_ts_.statevars()) {
      cache[to_prover_solver_.transfer_term(v)] = v;
    }
    for (auto v : orig_ts_.inputvars()) {
      cache[to_prover_solver_.transfer_term(v)] = v;
    }

    transfer_to_prover_as = [this](const Term & t, SortKind sk) {
      return to_prover_solver_.transfer_term(t, sk);
    };
    transfer_to_orig_ts_as = [&to_orig_ts_solver](const Term & t, SortKind sk) {
      return to_orig_ts_solver.transfer_term(t, sk);
    };
  }

  for (auto wit_map : witness_) {
    out.push_back(UnorderedTermMap());
    UnorderedTermMap & map = out.back();

    for (auto v : orig_ts_.statevars()) {
      SortKind sk = v->get_sort()->get_sort_kind();
      Term pv = transfer_to_prover_as(v, sk);
      map[v] = transfer_to_orig_ts_as(wit_map[pv], sk);
    }

    for (auto v : orig_ts_.inputvars()) {
      SortKind sk = v->get_sort()->get_sort_kind();
      Term pv = transfer_to_prover_as(v, sk);
      map[v] = transfer_to_orig_ts_as(wit_map[pv], sk);
    }

    for (auto elem : orig_ts_.named_terms()) {
      SortKind sk = elem.second->get_sort()->get_sort_kind();
      Term pt = transfer_to_prover_as(elem.second, sk);
      map[elem.second] = transfer_to_orig_ts_as(wit_map[pt], sk);
    }
  }

  return true;
}

Term Prover::invar()
{
  if (!invar_)
  {
    throw PonoException("Failed to return invar. Be sure that the property was proven "
                        "by an engine the supports returning invariants.");
  }
  return to_orig_ts(invar_, BOOL);
}

Term Prover::to_orig_ts(Term t, SortKind sk)
{
  if (solver_ == orig_ts_.solver()) {
    // don't need to transfer terms if the solvers are the same
    return t;
  } else {
    /* TODO: double-check that transferring terms still works as
       intended in this branch when COI is used. */
    if (options_.static_coi_)
      throw PonoException(
          "Temporary restriction: cone-of-influence analysis "
          "currently incompatible with witness generation.");
    // need to add symbols to cache
    TermTranslator to_orig_ts_solver(orig_ts_.solver());
    UnorderedTermMap & cache = to_orig_ts_solver.get_cache();
    for (auto v : orig_ts_.statevars()) {
      cache[to_prover_solver_.transfer_term(v)] = v;
      Term nv = orig_ts_.next(v);
      cache[to_prover_solver_.transfer_term(nv)] = v;
    }
    for (auto v : orig_ts_.inputvars()) {
      cache[to_prover_solver_.transfer_term(v)] = v;
    }
    return to_orig_ts_solver.transfer_term(t, sk);
  }
}

Term Prover::to_orig_ts(Term t)
{
  return to_orig_ts(t, t->get_sort()->get_sort_kind());
}

bool Prover::compute_witness()
{
  // TODO: make sure the solver state is SAT

  for (int i = 0; i <= reached_k_; ++i) {
    witness_.push_back(UnorderedTermMap());
    UnorderedTermMap & map = witness_.back();

    for (auto v : ts_.statevars()) {
      Term vi = unroller_.at_time(v, i);
      Term r = solver_->get_value(vi);
      map[v] = r;
    }

    for (auto v : ts_.inputvars()) {
      Term vi = unroller_.at_time(v, i);
      Term r = solver_->get_value(vi);
      map[v] = r;
    }

    for (auto elem : ts_.named_terms()) {
      Term ti = unroller_.at_time(elem.second, i);
      map[elem.second] = solver_->get_value(ti);
    }
  }

  return true;
}

}  // namespace pono
