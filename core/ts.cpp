/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan
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

#include "core/ts.h"

#include <functional>

#include "assert.h"
#include "smt-switch/substitution_walker.h"
#include "smt-switch/utils.h"

using namespace smt;
using namespace std;

namespace pono {

void swap(TransitionSystem & ts1, TransitionSystem & ts2)
{
  std::swap(ts1.solver_, ts2.solver_);
  std::swap(ts1.init_, ts2.init_);
  std::swap(ts1.trans_, ts2.trans_);
  std::swap(ts1.statevars_, ts2.statevars_);
  std::swap(ts1.next_statevars_, ts2.next_statevars_);
  std::swap(ts1.inputvars_, ts2.inputvars_);
  std::swap(ts1.named_terms_, ts2.named_terms_);
  std::swap(ts1.term_to_name_, ts2.term_to_name_);
  std::swap(ts1.state_updates_, ts2.state_updates_);
  std::swap(ts1.next_map_, ts2.next_map_);
  std::swap(ts1.curr_map_, ts2.curr_map_);
  std::swap(ts1.functional_, ts2.functional_);
  std::swap(ts1.deterministic_, ts2.deterministic_);
  std::swap(ts1.constraints_, ts2.constraints_);
}

TransitionSystem & TransitionSystem::operator=(TransitionSystem other)
{
  swap(*this, other);
  return *this;
}

TransitionSystem::TransitionSystem(const TransitionSystem & other_ts,
                                   TermTranslator & tt)
{
  function<Term(const Term &)> transfer;
  function<Term(const Term &, SortKind)> transfer_as;
  if (other_ts.solver() == tt.get_solver()) {
    // if the solvers are the same, don't need to transfer
    transfer = [](const Term & t) { return t; };
    // assume you don't need to do sort-casting for terms from the same solver
    transfer_as = [](const Term & t, SortKind sk) { return t; };
  } else {
    transfer = [&tt](const Term & t) { return tt.transfer_term(t); };
    transfer_as = [&tt](const Term & t, SortKind sk) {
      return tt.transfer_term(t, sk);
    };
  }

  solver_ = tt.get_solver();
  // transfer init and trans -- expect them to be boolean
  // will cast if underlying solver aliases Bool/BV1
  init_ = transfer_as(other_ts.init_, BOOL);
  trans_ = transfer_as(other_ts.trans_, BOOL);

  // populate data structures with translated terms

  for (const auto & v : other_ts.statevars_) {
    statevars_.insert(transfer(v));
  }

  for (const auto & v : other_ts.inputvars_) {
    inputvars_.insert(transfer(v));
  }

  for (const auto & v : other_ts.next_statevars_) {
    next_statevars_.insert(transfer(v));
  }

  for (const auto & elem : other_ts.named_terms_) {
    named_terms_[elem.first] = transfer(elem.second);
  }

  for (const auto & elem : other_ts.term_to_name_) {
    term_to_name_[transfer(elem.first)] = elem.second;
  }

  // variables might have already be in the TermTranslator cache
  // with a different sort (due to sort aliasing)
  // use the SortKind as a hint when transferring
  // sorts of the two terms should match for state updates and next_map
  Term key, val;
  for (const auto & elem : other_ts.state_updates_) {
    key = transfer(elem.first);
    val = transfer_as(elem.second, key->get_sort()->get_sort_kind());
    assert(key->get_sort() == val->get_sort());
    state_updates_[key] = val;
  }
  for (const auto & elem : other_ts.next_map_) {
    key = transfer(elem.first);
    val = transfer_as(elem.second, key->get_sort()->get_sort_kind());
    next_map_[key] = val;
  }

  for (const auto & elem : other_ts.curr_map_) {
    curr_map_[transfer(elem.first)] = transfer(elem.second);
  }

  /* Constraints collected in vector 'constraints_' were part of init_
     and/or trans_ and were transferred already above. Hence these
     terms should be in the term translator cache. */
  for (const auto & e : other_ts.constraints_) {
    constraints_.push_back({ transfer_as(e.first, BOOL), e.second });
  }
  functional_ = other_ts.functional_;
  deterministic_ = other_ts.deterministic_;
}

bool TransitionSystem::operator==(const TransitionSystem & other) const
{
  return (solver_ == other.solver_ &&
          init_ == other.init_ &&
          trans_ == other.trans_ &&
          statevars_ == other.statevars_ &&
          next_statevars_ == other.next_statevars_ &&
          inputvars_ == other.inputvars_ &&
          named_terms_ == other.named_terms_ &&
          term_to_name_ == other.term_to_name_ &&
          state_updates_ == other.state_updates_ &&
          next_map_ == other.next_map_ &&
          curr_map_ == other.curr_map_ &&
          functional_ == other.functional_ &&
          deterministic_ == other.deterministic_ &&
          constraints_ == other.constraints_);
}

bool TransitionSystem::operator!=(const TransitionSystem & other) const
{
  return !(*this == other);
}

void TransitionSystem::set_init(const Term & init)
{
  // TODO: only do this check in debug mode
  if (!only_curr(init)) {
    throw PonoException(
        "Initial state constraints should only use current state variables");
  }

  init_ = init;
}

void TransitionSystem::constrain_init(const Term & constraint)
{
  // TODO: Only do this check in debug mode
  if (!only_curr(constraint)) {
    throw PonoException(
        "Initial state constraints should only use current state variables");
  }
  init_ = solver_->make_term(And, init_, constraint);
}

void TransitionSystem::assign_next(const Term & state, const Term & val)
{
  // TODO: only do this check in debug mode
  if (statevars_.find(state) == statevars_.end()) {
    throw PonoException("Unknown state variable");
  }

  if (!no_next(val)) {
    throw PonoException(
        "Got a symbolic that is not a current state or input variable in RHS "
        "of functional assignment");
  }

  if (state_updates_.find(state) != state_updates_.end()) {
    throw PonoException("State variable " + state->to_string()
                        + " already has next-state logic assigned.");
  }

  state_updates_[state] = val;
  trans_ = solver_->make_term(
      And, trans_, solver_->make_term(Equal, next_map_.at(state), val));

  // if not functional, then we cannot guarantee deterministm
  // if it is functional, depends on if all state variables
  // have updates
  // technically not even functional if there are constraints
  // TODO: revisit this and possibly rename functional/deterministic
  if (functional_ && !constraints_.size()) {
    deterministic_ = (state_updates_.size() == statevars_.size());
  }
}

void TransitionSystem::add_invar(const Term & constraint)
{
  // invariants can make it so not every state has a next state
  // TODO: revisit this and possibly rename functional/deterministic
  deterministic_ = false;

  // TODO: only check this in debug mode
  if (only_curr(constraint)) {
    init_ = solver_->make_term(And, init_, constraint);
    trans_ = solver_->make_term(And, trans_, constraint);
    Term next_constraint = solver_->substitute(constraint, next_map_);
    // add the next-state version
    trans_ = solver_->make_term(And, trans_, next_constraint);
    constraints_.push_back({ constraint, true });
  } else {
    throw PonoException("Invariants should be over current states only.");
  }
}

void TransitionSystem::constrain_inputs(const Term & constraint)
{
  // constraints can make it so not every state has a next state
  // TODO: revisit this and possibly rename functional/deterministic
  deterministic_ = false;

  if (no_next(constraint)) {
    trans_ = solver_->make_term(And, trans_, constraint);
    constraints_.push_back({ constraint, true });
  } else {
    throw PonoException("Cannot have next-states in an input constraint.");
  }
}

void TransitionSystem::add_constraint(const Term & constraint,
                                      bool to_init_and_next)
{
  // constraints can make it so not every state has a next state
  // TODO: revisit this and possibly rename functional/deterministic
  deterministic_ = false;

  if (only_curr(constraint)) {
    trans_ = solver_->make_term(And, trans_, constraint);

    if (to_init_and_next) {
      init_ = solver_->make_term(And, init_, constraint);
      Term next_constraint = solver_->substitute(constraint, next_map_);
      trans_ = solver_->make_term(And, trans_, next_constraint);
    }
    constraints_.push_back({ constraint, to_init_and_next });
  } else if (no_next(constraint)) {
    trans_ = solver_->make_term(And, trans_, constraint);
    constraints_.push_back({ constraint, to_init_and_next });
  } else {
    throw PonoException("Constraint cannot have next states");
  }
}

void TransitionSystem::name_term(const string name, const Term & t)
{
  auto it = named_terms_.find(name);
  if (it != named_terms_.end() && t != it->second) {
    throw PonoException("Name " + name + " has already been used.");
  }
  named_terms_[name] = t;
  // save this name as a representative (might overwrite)
  term_to_name_[t] = name;
}

Term TransitionSystem::make_inputvar(const string name, const Sort & sort)
{
  Term input = solver_->make_symbol(name, sort);
  add_inputvar(input);
  return input;
}

Term TransitionSystem::make_statevar(const string name, const Sort & sort)
{
  // set to false until there is a next state update for this statevar
  deterministic_ = false;

  Term state = solver_->make_symbol(name, sort);
  Term next_state = solver_->make_symbol(name + ".next", sort);
  add_statevar(state, next_state);
  return state;
}

Term TransitionSystem::curr(const Term & term) const
{
  return solver_->substitute(term, curr_map_);
}

Term TransitionSystem::next(const Term & term) const
{
  if (next_map_.find(term) != next_map_.end()) {
    return next_map_.at(term);
  }
  return solver_->substitute(term, next_map_);
}

bool TransitionSystem::is_curr_var(const Term & sv) const
{
  return (statevars_.find(sv) != statevars_.end());
}

bool TransitionSystem::is_next_var(const Term & sv) const
{
  return (next_statevars_.find(sv) != next_statevars_.end());
}

bool TransitionSystem::is_input_var(const Term & sv) const
{
  return (inputvars_.find(sv) != inputvars_.end());
}

std::string TransitionSystem::get_name(const Term & t) const
{
  const auto & it = term_to_name_.find(t);
  if (it != term_to_name_.end()) {
    return it->second;
  }
  return t->to_string();
}

smt::Term TransitionSystem::lookup(std::string name) const
{
  const auto & it = named_terms_.find(name);
  if (it == named_terms_.end()) {
    throw PonoException("Could not find term named: " + name);
  }
  return it->second;
}

void TransitionSystem::add_statevar(const Term & cv, const Term & nv)
{
  // TODO: this runs even if called from make_statevar
  //       could refactor entirely, or just pass a flag
  //       saying whether to check these things or not

  if (statevars_.find(cv) != statevars_.end()) {
    throw PonoException("Cannot redeclare a state variable");
  }

  if (next_statevars_.find(nv) != next_statevars_.end()) {
    throw PonoException("Cannot redeclare a state variable");
  }

  if (next_statevars_.find(cv) != next_statevars_.end()) {
    throw PonoException(
        "Cannot use an existing next state variable as a current state var");
  }

  if (statevars_.find(nv) != statevars_.end()) {
    throw PonoException(
        "Cannot use an existing state variable as a next state var");
  }

  // if using an input variable, remove from set
  // will be a state variable now
  if (inputvars_.find(cv) != inputvars_.end()) {
    bool success = inputvars_.erase(cv);
    assert(success);
  }

  if (inputvars_.find(nv) != inputvars_.end()) {
    bool success = inputvars_.erase(nv);
    assert(success);
  }

  statevars_.insert(cv);
  next_statevars_.insert(nv);
  next_map_[cv] = nv;
  curr_map_[nv] = cv;
  // automatically include in named_terms
  name_term(cv->to_string(), cv);
  name_term(nv->to_string(), nv);
}

void TransitionSystem::add_inputvar(const Term & v)
{
  // TODO: this check is running even when used by make_inputvar
  //       could refactor entirely or just pass a boolean saying whether or not
  //       to check these things
  if (statevars_.find(v) != statevars_.end()
      || next_statevars_.find(v) != next_statevars_.end()
      || inputvars_.find(v) != inputvars_.end()) {
    throw PonoException(
        "Cannot reuse an existing variable as an input variable");
  }

  inputvars_.insert(v);
  // automatically include in named_terms
  name_term(v->to_string(), v);
}

// term building methods -- forwards to SmtSolver solver_

Sort TransitionSystem::make_sort(const std::string name, uint64_t arity)
{
  return solver_->make_sort(name, arity);
}

Sort TransitionSystem::make_sort(const SortKind sk)
{
  return solver_->make_sort(sk);
}

Sort TransitionSystem::make_sort(const SortKind sk, uint64_t size)
{
  return solver_->make_sort(sk, size);
}

Sort TransitionSystem::make_sort(const SortKind sk, const Sort & sort1)
{
  return solver_->make_sort(sk, sort1);
}

Sort TransitionSystem::make_sort(const SortKind sk,
                                 const Sort & sort1,
                                 const Sort & sort2)
{
  return solver_->make_sort(sk, sort1, sort2);
}

Sort TransitionSystem::make_sort(const SortKind sk,
                                 const Sort & sort1,
                                 const Sort & sort2,
                                 const Sort & sort3)
{
  return solver_->make_sort(sk, sort1, sort2, sort3);
}

Sort TransitionSystem::make_sort(const SortKind sk, const SortVec & sorts)
{
  return solver_->make_sort(sk, sorts);
}

Term TransitionSystem::make_term(bool b) { return solver_->make_term(b); }

Term TransitionSystem::make_term(int64_t i, const Sort & sort)
{
  return solver_->make_term(i, sort);
}

Term TransitionSystem::make_term(const std::string val,
                                 const Sort & sort,
                                 uint64_t base)
{
  return solver_->make_term(val, sort, base);
}

Term TransitionSystem::make_term(const Term & val, const Sort & sort)
{
  return solver_->make_term(val, sort);
}

Term TransitionSystem::make_term(const Op op, const Term & t)
{
  return solver_->make_term(op, t);
}

Term TransitionSystem::make_term(const Op op, const Term & t0, const Term & t1)
{
  return solver_->make_term(op, t0, t1);
}

Term TransitionSystem::make_term(const Op op,
                                 const Term & t0,
                                 const Term & t1,
                                 const Term & t2)
{
  return solver_->make_term(op, t0, t1, t2);
}

Term TransitionSystem::make_term(const Op op, const TermVec & terms)
{
  return solver_->make_term(op, terms);
}

void TransitionSystem::rebuild_trans_based_on_coi(
    const UnorderedTermSet & state_vars_in_coi,
    const UnorderedTermSet & input_vars_in_coi)
{
  /* Clear current transition relation 'trans_'. */
  trans_ = solver_->make_term(true);

  /* Add next-state functions for state variables in COI. */
  for (const auto & state_var : state_vars_in_coi) {
    Term next_func = NULL;
    const auto & elem = state_updates_.find(state_var);
    if (elem != state_updates_.end())
      next_func = elem->second;
    /* May find state variables without next-function. */
    if (next_func != NULL) {
        Term eq = solver_->make_term(Equal, next_map_.at(state_var), next_func);
        trans_ = solver_->make_term(And, trans_, eq);
      }
  }

  /* Add global constraints added to previous 'trans_'. */
  // TODO: check potential optimizations in removing global constraints
  std::vector<std::pair<smt::Term, bool>> prev_constraints = constraints_;
  constraints_.clear();
  for (const auto & e : prev_constraints) {
    add_constraint(e.first, e.second);
  }

  statevars_.clear();
  // Have to add any state variables in init back in
  // this is because COI doesn't consider init and if
  // we remove those state variables then the TS is
  // ill-formed (e.g. init will contain unknown symbols)
  // this shouldn't affect performance much, because
  // variables in init that are not in the COI *only*
  // appear in init
  get_free_symbolic_consts(init_, statevars_);
  for (const auto & var : state_vars_in_coi) {
    statevars_.insert(var);
  }

  inputvars_.clear();
  for (const auto & var : input_vars_in_coi) {
    inputvars_.insert(var);
  }

  smt::UnorderedTermMap reduced_state_updates;
  for (const auto & var : state_vars_in_coi) {
    const auto & elem = state_updates_.find(var);
    if (elem != state_updates_.end()) {
      Term next_func = elem->second;
      reduced_state_updates[var] = next_func;
    }
  }
  state_updates_ = reduced_state_updates;

  /* update named_terms and term_to_name_ by removing terms that
     no longer exist in the system
   */
  unordered_map<string, Term> reduced_named_terms;
  unordered_map<Term, string> reduced_term_to_name;
  UnorderedTermSet free_vars;
  for (const auto & elem : named_terms_) {
    free_vars.clear();
    get_free_symbolic_consts(elem.second, free_vars);
    bool all_in_sys = true;
    Term currvar;
    for (const auto & v : free_vars) {
      // v is an input variable, current variable, or next variable
      // we want the current version of a state variable
      const auto & it = curr_map_.find(v);
      if (it != curr_map_.end()) {
        // get the current state version of a next variable
        currvar = it->second;
      } else {
        currvar = v;
      }

      if (statevars_.find(currvar) == statevars_.end()
          && inputvars_.find(currvar) == inputvars_.end()) {
        all_in_sys = false;
        break;
      }
    }

    if (all_in_sys) {
      reduced_named_terms[elem.first] = elem.second;
      // NOTE: name might not be the same as elem.first
      //       need to use the representative name
      //       stored in term_to_name_
      reduced_term_to_name[elem.second] = term_to_name_.at(elem.second);
    }
  }
  named_terms_ = reduced_named_terms;
  term_to_name_ = reduced_term_to_name;
}

// protected methods

bool TransitionSystem::contains(const Term & term,
                                UnorderedTermSetPtrVec term_sets) const
{
  UnorderedTermSet visited;
  TermVec to_visit{ term };
  Term t;
  while (to_visit.size()) {
    t = to_visit.back();
    to_visit.pop_back();

    if (visited.find(t) != visited.end()) {
      // cache hit
      continue;
    }

    if (t->is_symbolic_const()) {
      bool in_atleast_one = false;
      for (const auto & ts : term_sets) {
        if (ts->find(t) != ts->end()) {
          in_atleast_one = true;
          break;
        }
      }

      if (!in_atleast_one) {
        return false;
      }
    }

    visited.insert(t);
    for (const auto & c : t) {
      to_visit.push_back(c);
    }
  }

  return true;
}

bool TransitionSystem::only_curr(const Term & term) const
{
  return contains(term, UnorderedTermSetPtrVec{ &statevars_ });
}

bool TransitionSystem::no_next(const Term & term) const
{
  return contains(term, UnorderedTermSetPtrVec{ &statevars_, &inputvars_ });
}

void TransitionSystem::drop_state_updates(const TermVec & svs)
{
  for (const auto & sv : svs) {
    if (!is_curr_var(sv)) {
      throw PonoException("Got non-state var in drop_state_updates");
    }
    state_updates_.erase(sv);
  }

  // now rebuild trans
  /* Clear current transition relation 'trans_'. */
  trans_ = solver_->make_term(true);

  /* Add next-state functions for state variables in COI. */
  for (const auto & elem : state_updates_) {
    assert(elem.second);  // should be non-null if in map
    Term eq = solver_->make_term(Equal, next_map_.at(elem.first), elem.second);
    trans_ = solver_->make_term(And, trans_, eq);
  }

  /* Add global constraints added to previous 'trans_'. */
  for (const auto & e : constraints_) {
    add_constraint(e.first, e.second);
  }
}

void TransitionSystem::promote_inputvar(const Term & iv)
{
  size_t num_erased = inputvars_.erase(iv);
  if (!num_erased) {
    throw PonoException("Tried to promote non-input to state variable: "
                        + iv->to_string());
  }

  // now turn into a state variable

  // set to false until there is a next state update for this statevar
  deterministic_ = false;
  Term next_state_var =
      solver_->make_symbol(iv->to_string() + ".next", iv->get_sort());
  add_statevar(iv, next_state_var);
}

void TransitionSystem::replace_terms(const UnorderedTermMap & to_replace)
{
  // first check that all the replacements contain known symbols
  UnorderedTermSetPtrVec all_symbols(
      { &statevars_, &inputvars_, &next_statevars_ });
  for (const auto & elem : to_replace) {
    bool known = contains(elem.first, all_symbols);
    known &= contains(elem.second, all_symbols);
    if (!known) {
      throw PonoException("Got an unknown symbol in replace_terms map");
    }
  }

  // use a substitution walker because
  //   1. it keeps a persistent cache
  //   2. it supports substituting arbitrary terms (e.g. not just mapping from
  //   symbols)
  SubstitutionWalker sw(solver_, to_replace);

  // now rebuild terms in every data structure with replacements
  init_ = sw.visit(init_);
  if (!only_curr(init_)) {
    throw PonoException(
        "Replaced a state variable appearing in init with an input in "
        "replace_terms");
  }
  trans_ = sw.visit(trans_);

  unordered_map<string, Term> new_named_terms;
  unordered_map<Term, string> new_term_to_name;
  for (auto elem : named_terms_) {
    new_named_terms[elem.first] = sw.visit(elem.second);
    new_term_to_name[sw.visit(elem.second)] = term_to_name_.at(elem.second);
  }
  named_terms_ = new_named_terms;
  term_to_name_ = new_term_to_name;

  // NOTE: don't need to update vars, let COI reduction handle that
  UnorderedTermMap new_state_updates;
  Term sv, update;
  for (auto elem : state_updates_) {
    sv = elem.first;
    sv = sw.visit(sv);
    update = sw.visit(elem.second);
    if (functional_ && !no_next(update)) {
      throw PonoException(
          "Got a next state variable in a state update for a functional "
          "TransitionSystem in replace_terms");
    }
    new_state_updates[sv] = update;
  }
  state_updates_ = new_state_updates;

  UnorderedTermMap new_next_map_;
  UnorderedTermMap new_curr_map_;
  Term c, n;
  for (const auto & elem : next_map_) {
    c = elem.first;
    n = elem.second;
    c = sw.visit(c);
    n = sw.visit(n);
    new_next_map_[c] = n;
    assert(curr_map_.at(elem.second) == elem.first);
    new_curr_map_[n] = c;
  }
  next_map_ = new_next_map_;
  curr_map_ = new_curr_map_;

  vector<pair<Term, bool>> new_constraints;
  new_constraints.reserve(constraints_.size());
  for (const auto & e : constraints_) {
    new_constraints.push_back(e);
  }
  constraints_ = new_constraints;
}

bool TransitionSystem::known_symbols(const Term & term) const
{
  return contains(
      term,
      UnorderedTermSetPtrVec{ &statevars_, &inputvars_, &next_statevars_ });
}

}  // namespace pono
