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
    compute_coi();
    ts_.rebuild_trans_based_on_coi(statevars_in_coi_);
  }
}

ProverResult Prover::prove() { return check_until(INT_MAX); }

//TODO: mainly for debugging but could make it an optional output at
//certain verbosity level, check logger usage
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

//TODO: mainly for debugging but could make it an optional output at
//certain verbosity level, check logger usage
void Prover::print_bad_property_coi()
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

}

void Prover::collect_coi_term(UnorderedTermSet & set, const Term & term)
{
  if (set.find(term) == set.end())
    set.insert(term);
}
  
void Prover::compute_term_coi(const Term & term,
                              UnorderedTermSet & new_coi_state_vars,
                              UnorderedTermSet & new_coi_input_vars)
{
  assert(term != NULL);
  TermVec open_terms;
  open_terms.push_back (term);
  Term cur;

  /* Depth-first search of term structure 'term'. Collect all
     encountered state and input variables and store them in
     'new_coi_state_vars' and 'new_coi_input_vars'. */
  while (!open_terms.empty()) {
    cur = open_terms.back();
    open_terms.pop_back ();

    /* Use global set 'coi_visited_terms' here to avoid visiting terms
       multiple times when we call this function on different terms. */
    if (coi_visited_terms_.find(cur) == coi_visited_terms_.end()) {
      /* Cache 'cur' and push its children. */
      coi_visited_terms_.insert(cur);
      cout << "  visiting COI term: " << cur << "\n";
      //TODO: check in smt-switch if '->is_symbol' captures statevars properly
      if (cur->is_symbol()) {
          cout << "    ..is symbol\n";
          if (ts_.statevars().find(cur) != ts_.statevars().end()) {
            assert(ts_.inputvars().find(cur) == ts_.inputvars().end());
            assert(!ts_.is_next_var(cur));
            cout << "collect COI statevar " << cur << "\n";
            collect_coi_term(new_coi_state_vars, cur);
          }
          else if (ts_.inputvars().find(cur) != ts_.inputvars().end()) {
            assert(!ts_.is_curr_var(cur));
            assert(!ts_.is_next_var(cur));
            cout << "collect COI inputvar " << cur << "\n";
            collect_coi_term(new_coi_input_vars, cur);
          }
        }
      
      for (auto child : cur) {
        //        cout << "    pushing child: " << child << "\n";
        open_terms.push_back(child);
      }
    }
  }
}

void Prover::compute_coi_next_state_funcs(UnorderedTermSet & new_coi_state_vars,
                                          UnorderedTermSet & new_coi_input_vars)
{  
  cout << "COI analysis: next-state functions" << endl;

  TermVec unprocessed_state_vars;
  for (auto state_var : statevars_in_coi_)
    unprocessed_state_vars.push_back(state_var);

  while (!unprocessed_state_vars.empty()) {
    Term state_var = unprocessed_state_vars.back();
    unprocessed_state_vars.pop_back();
    assert(ts_.is_curr_var(state_var));
    const smt::UnorderedTermMap state_updates = ts_.state_updates();
    Term next_func = NULL;
    /* May find state variables without next-function. */
    auto elem = ts_.state_updates().find(state_var);
    if (elem != ts_.state_updates().end())
      next_func = elem->second;

    assert(new_coi_state_vars.empty());
    assert(new_coi_input_vars.empty());

    if (next_func != NULL)
      compute_term_coi(next_func, new_coi_state_vars, new_coi_input_vars);

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
  
void Prover::compute_coi()
{
  assert (coi_visited_terms_.empty());
  //TODO: remove printing or toggle wrt verbosity level
  //print_bad_property_coi();

  cout << "Starting COI analysis" << endl;

  UnorderedTermSet new_coi_state_vars;
  UnorderedTermSet new_coi_input_vars;

  cout << "Starting COI analysis: bad-term" << endl;
  compute_term_coi(bad_, new_coi_state_vars, new_coi_input_vars);

  assert (statevars_in_coi_.empty());
  assert (inputvars_in_coi_.empty());

  for (auto sv : new_coi_state_vars)
    statevars_in_coi_.insert(sv);
  for (auto sv : new_coi_input_vars)
    inputvars_in_coi_.insert(sv);

  new_coi_state_vars.clear();
  new_coi_input_vars.clear();  

  compute_coi_next_state_funcs(new_coi_state_vars, new_coi_input_vars);

  //TODO/NOTE: we do NOT traverse 'init_' constraint to search for new
  //state vars; we leave 'init_' as is with COI analysis and only
  //rebuild 'trans_'; 'init_' can be any constraint and could be
  //difficult to figure out which parts are in COI; also, any state
  //vars in 'init_' that have not been identified as part of the COI
  //by checking 'bad_' and the next-state functions cannot have an
  //influence on the property 'bad_'.
  
  cout << "COI analysis completed" << "\n";
  for (auto var : statevars_in_coi_)
    cout << "  found COI statevar " << var << "\n";
  for (auto var : inputvars_in_coi_)
    cout << "  found COI inputvar " << var << "\n";

  cout << "Original system had:" << "\n";
  for (auto var : ts_.statevars())
    cout << "  statevar " << var << "\n";
  for (auto var : ts_.inputvars())
    cout << "  inputvar " << var << "\n";
}
  
bool Prover::witness(std::vector<UnorderedTermMap> & out)
{
  // TODO: make sure the solver state is SAT

  function<Term(const Term &, SortKind)> transfer_to_prover_as;
  function<Term(const Term &, SortKind)> transfer_to_orig_ts_as;
  TermTranslator to_orig_ts_solver(orig_ts_.solver());
  if (solver_ == orig_ts_.solver()) {
    // don't need to transfer terms if the solvers are the same
    transfer_to_prover_as = [](const Term & t, SortKind sk) { return t; };
    transfer_to_orig_ts_as = [](const Term & t, SortKind sk) { return t; };
  } else {
    /* TODO: double-check that transferring terms still works as
       intended when COI is used. */
    if (options_.static_coi_)
      throw PonoException("Temporary restriction: cone-of-influence analysis "\
                          "currently incompatible with transferring terms.");
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

  for (int i = 0; i <= reached_k_; ++i) {
    out.push_back(UnorderedTermMap());
    UnorderedTermMap & map = out.back();

    for (auto v : orig_ts_.statevars()) {
      SortKind sk = v->get_sort()->get_sort_kind();
      Term vi = unroller_.at_time(transfer_to_prover_as(v, sk), i);
      Term r = solver_->get_value(vi);
      map[v] = transfer_to_orig_ts_as(r, sk);
    }

    for (auto v : orig_ts_.inputvars()) {
      SortKind sk = v->get_sort()->get_sort_kind();
      Term vi = unroller_.at_time(transfer_to_prover_as(v, sk), i);
      Term r = solver_->get_value(vi);
      map[v] = transfer_to_orig_ts_as(r, sk);
    }

    for (auto elem : orig_ts_.named_terms()) {
      SortKind sk = elem.second->get_sort()->get_sort_kind();
      Term ti = unroller_.at_time(transfer_to_prover_as(elem.second, sk), i);
      map[elem.second] = transfer_to_orig_ts_as(solver_->get_value(ti), sk);
    }
  }

  return true;
}

}  // namespace pono
