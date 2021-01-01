/*********************                                                  */
/*! \file ic3sa.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief IC3 with Syntax-Guided Abstraction based on
**
**        Model Checking of Verilog RTL Using IC3 with Syntax-Guided
*Abstraction.
**            -- Aman Goel, Karem Sakallah
**
**
**  within Pono, we are building on the bit-level IC3 instead of directly
**  on IC3Base, because a lot of the functionality is the same
**  In particular, we don't need to override inductive generalization
**
**/

#include "engines/ic3sa.h"

#include "assert.h"
#include "smt-switch/utils.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"
#include "utils/term_walkers.h"

using namespace smt;
using namespace std;

// TODO add option to use interpolants in inductive generalization

namespace pono {

// main IC3SA implementation

IC3SA::IC3SA(Property & p, const smt::SmtSolver & solver, PonoOptions opt)
    : super(p, solver, opt), fcoi_(*ts_, 0)
{
}

IC3Formula IC3SA::get_model_ic3formula(TermVec * out_inputs,
                                       TermVec * out_nexts) const
{
  TermVec cube_lits;
  // first populate with predicates
  for (const auto & p : predset_) {
    if (solver_->get_value(p) == solver_true_) {
      cube_lits.push_back(p);
    } else {
      cube_lits.push_back(solver_->make_term(Not, p));
    }
  }

  // TODO make sure that projecting on state variables here makes sense
  EquivalenceClasses ec = get_equivalence_classes_from_model(ts_->statevars());
  construct_partition(ec, cube_lits);
  IC3Formula cube = ic3formula_conjunction(cube_lits);
  assert(ic3formula_check_valid(cube));

  if (out_nexts) {
    for (const auto & sv : ts_->statevars()) {
      out_nexts->push_back(
          solver_->make_term(Equal, sv, solver_->get_value(sv)));
    }
  }

  if (out_inputs) {
    for (const auto & iv : ts_->inputvars()) {
      out_inputs->push_back(
          solver_->make_term(Equal, iv, solver_->get_value(iv)));
    }
  }

  return cube;
}

bool IC3SA::ic3formula_check_valid(const IC3Formula & u) const
{
  Sort boolsort = solver_->make_sort(BOOL);
  // check that children are literals
  Op op;
  for (auto c : u.children) {
    op = c->get_op();

    if (c->get_sort() != boolsort) {
      return false;
    }

    // include bvnot for boolector
    if (op == Not || op == BVNot) {
      c = smart_not(c);
    }

    // an equality is a predicate, so we just need to check
    // is_predicate, not specifically for equalities
    if (!is_predicate(c, boolsort, true)) {
      return false;
    }
  }

  // for now not checking the actual term (e.g. u.term)
  // it's somewhat hard if the underlying solver does rewriting

  // got through all checks without failing
  return true;
}

IC3Formula IC3SA::generalize_predecessor(size_t i, const IC3Formula & c)
{
  // TODO: use the JustifyCOI algorithm from the paper
  //       e.g. partial_model
  // compute cone-of-influence of target c
  fcoi_.compute_coi({ c.term });
  const UnorderedTermSet & coi_symbols = fcoi_.statevars_in_coi();
  assert(coi_symbols.size() <= ts_->statevars().size());

  logger.log(
      2,
      "IC3SA::generalize_predecessor projecting on {}/{} state variables",
      coi_symbols.size(),
      ts_->statevars().size());

  TermVec cube_lits;

  // first populate with predicates
  UnorderedTermSet free_symbols;
  for (const auto & p : predset_) {
    free_symbols.clear();
    get_free_symbolic_consts(p, free_symbols);
    for (const auto & fv : free_symbols) {
      if (coi_symbols.find(fv) != coi_symbols.end()) {
        continue;
      }
    }

    if (solver_->get_value(p) == solver_true_) {
      cube_lits.push_back(p);
    } else {
      cube_lits.push_back(solver_->make_term(Not, p));
    }
  }

  EquivalenceClasses ec = get_equivalence_classes_from_model(coi_symbols);
  construct_partition(ec, cube_lits);
  IC3Formula cube = ic3formula_conjunction(cube_lits);
  assert(ic3formula_check_valid(cube));
  return cube;
}

void IC3SA::check_ts() const
{
  // TODO: everything in ic3sa should be a state variable (at least according to
  // paper)
  //       might work if we just remove input variables from the subterms
  // TODO: add option to promote all inputs to be state vars
  // TODO: add support for arrays

  for (const auto & sv : ts_->statevars())
  {
    SortKind sk = sv->get_sort()->get_sort_kind();
    if (sk != BOOL && sk != BV)
    {
      throw PonoException("IC3SA currently only supports bit-vectors");
    }
  }
  for (const auto & iv : ts_->inputvars())
  {
    SortKind sk = iv->get_sort()->get_sort_kind();
    if (sk != BOOL && sk != BV)
    {
      throw PonoException("IC3SA currently only supports bit-vectors");
    }
  }
}

RefineResult IC3SA::refine()
{
  // TODO try substituting in trans vs symbolic post-image
  //      trans substitution idea is hand-wavey, but want to end up
  //      with (cnt + 1 + 1) instead of just cnt = 1 for example
  // TODO add option for using interpolants
  assert(solver_context_ == 0);

  // recover the counterexample trace
  assert(check_intersects_initial(cex_pg_->target.term));
  vector<IC3Formula> cex;
  const ProofGoal * tmp = cex_pg_;
  while (tmp) {
    cex.push_back(tmp->target);
    assert(!tmp->target.disjunction);  // expecting a conjunction
    assert(ts_->only_curr(tmp->target.term));
    tmp = tmp->next;
  }

  size_t cex_length = cex.size();
  logger.log(1, "IC3SA::refine cex of length {}", cex_length);
  assert(cex_length);
  if (cex_length == 1) {
    // TODO make sure that this is indeed a concrete CEX
    return REFINE_NONE;
  }

  const UnorderedTermMap & state_updates = ts_->state_updates();

  // used to get rid of "unrolled" variables
  // after successfully refining
  UnorderedTermMap last_model_vals;

  // set up initial substitution map
  UnorderedTermMap subst;

  // assumps is for p_{i-1} /\ c_{i-1} /\ c_i' from paper
  TermVec assumps(cex[0].children.begin(), cex[0].children.end());
  assumps.push_back(ts_->init());
  Term trans = ts_->trans();
  Result r;
  bool refined = false;
  Term axiom;
  for (size_t i = 1; i < cex_length; ++i) {
    // add ci'
    // TODO use substitute_terms to save time when copying map
    for (const auto & c : cex[i].children) {
      assumps.push_back(ts_->next(c));
    }

    assert(solver_context_ == 0);
    push_solver_context();

    solver_->assert_formula(trans);

    for (const auto & a : assumps) {
      solver_->assert_formula(a);
    }

    r = solver_->check_sat();
    if (r.is_sat()) {
      assumps = symbolic_post_image(i, subst, last_model_vals);
      // add cube to it
      assumps.insert(
          assumps.end(), cex[i].children.begin(), cex[i].children.end());
    } else {
      refined = true;

      // get MUS
      // TODO check we're reducing over the correct vector
      //      should we include the cex cubes?
      TermVec unsatcore;
      reducer_.reduce_assump_unsatcore(trans, assumps, unsatcore);
      assert(unsatcore.size());
      // TODO consider removing ITEs at this step also
      Term m = make_and(unsatcore);

      // get rid of fresh symbolic constants for unconstrained variables
      // e.g. inputs and state variables with no state update
      // TODO consider trying to replace input variables by substitution
      //      I guess that means other terms that had the same value in
      //      the last satisfiable solver call? Not totally clear on that.
      m = solver_->substitute(m, last_model_vals);

      // instead of sv -> state_update, maps sv' -> state_update
      UnorderedTermMap next_updates;
      for (const auto & elem : state_updates) {
        next_updates[ts_->next(elem.first)] = elem.second;
      }
      // replace next-state variables with functional substitution
      m = solver_->substitute(m, next_updates);

      axiom = solver_->make_term(Not, m);
      assert(ts_->no_next(axiom));
      ts_->add_constraint(axiom);
      logger.log(2, "IC3SA::refine learning axiom: {}", axiom);

      // add all state variables to projection set
      UnorderedTermSet free_vars;
      get_free_symbolic_consts(axiom, free_vars);
      for (const auto & fv : free_vars) {
        if (ts_->is_curr_var(fv)) {
          projection_set_.insert(fv);
        }
      }

      bool new_terms_added = add_to_term_abstraction(axiom);
      assert(new_terms_added);
    }

    pop_solver_context();

    if (refined) {
      // expecting only one axioms
      // needs to be changed if we don't break from the loop upon finding an
      // axiom add semantics for trans_label_
      assert(axiom);
      assert(solver_context_ == 0);
      solver_->assert_formula(solver_->make_term(Implies, trans_label_, axiom));
      if (ts_->only_curr(axiom)) {
        solver_->assert_formula(
            solver_->make_term(Implies, trans_label_, ts_->next(axiom)));
      }
      break;
    }
  }

  // need to clear proof goals
  // cubes in proof goals aren't precise wrt to updated
  // terms (if refined)
  // so the partial trace it holds is no longer precise
  proof_goals_.clear();

  if (!refined) {
    // concrete counterexample
    return REFINE_NONE;
  }

  // TODO loop up to cex_length
  //      with implicit unrolling
  //      need to understand how T is being used
  //      with a functional unrolling in algorithm

  // TODO use symbolic_post_image
  // until the query becomes unsat
  // then add terms to term abstraction
  // (after substituting for inputs) and untiming
  // NOTE: seems easier to not use functional unroller
  //       need symbolic post-image *under current model*
  // TODO maybe have option for functional unroller
  // to not use @0 if never using other state variables
  // TODO figure out if we should project
  //      / how we limit the number of added terms
  // TODO get minimal unsat core
  // TODO add symbols from the MUS to the projection set permanently
  return REFINE_SUCCESS;
}

bool IC3SA::intersects_bad()
{
  push_solver_context();
  // assert the last frame (conjunction over clauses)
  assert_frame_labels(reached_k_ + 1);
  // see if it intersects with bad
  solver_->assert_formula(bad_);
  Result r = check_sat();

  if (r.is_sat()) {
    // start with a structural COI for reduction

    TermVec cube_lits;

    // first populate with predicates
    for (const auto & p : predset_) {
      if (!in_projection(p, vars_in_bad_)) {
        continue;
      }

      if (solver_->get_value(p) == solver_true_) {
        cube_lits.push_back(p);
      } else {
        cube_lits.push_back(solver_->make_term(Not, p));
      }
    }

    // TODO make sure that projecting on variables in bad here makes sense
    EquivalenceClasses ec = get_equivalence_classes_from_model(vars_in_bad_);
    construct_partition(ec, cube_lits);
    // reduce cube_lits
    TermVec red_c;
    reducer_.reduce_assump_unsatcore(smart_not(bad_), cube_lits, red_c);
    add_proof_goal(ic3formula_conjunction(red_c), reached_k_ + 1, NULL);
  }

  pop_solver_context();

  assert(!r.is_unknown());
  return r.is_sat();
}

void IC3SA::initialize()
{
  super::initialize();

  // IC3SA assumes input variables are modeled as state variables
  // with no update function
  // This seems important because otherwise we need to drop terms
  // containing input variables from IC3Formulas
  // because of destructive update, get copy of input variables first
  // TODO: decide whether to include this or not
  // UnorderedTermSet inputvars = ts_->inputvars();
  // for (const auto & iv : inputvars) {
  //   ts_->promote_inputvar(iv);
  // }
  // TODO with change above, don't need to check only_curr everywhere
  // technically should remove those checks (or at least guard with an assert)
  // assert(!ts_->inputvars().size());

  // set up initial term abstraction by getting all subterms
  // TODO consider starting with only a subset -- e.g. variables
  // TODO consider keeping a cache from terms to their free variables
  //      for use in COI
  SubTermCollector stc(solver_);
  stc.collect_subterms(ts_->init());
  stc.collect_subterms(ts_->trans());
  stc.collect_subterms(bad_);

  // TODO make sure this is right
  // I think we'll always project models onto at least state variables
  // so, we should prune those terms now
  // otherwise we'll do unnecessary iteration over them every time we get a
  // model

  for (const auto & elem : stc.get_subterms()) {
    const Sort & sort = elem.first;
    for (const auto & term : elem.second) {
      if (ts_->only_curr(term)) {
        term_abstraction_[sort].insert(term);
      }
    }
  }

  for (const auto & p : stc.get_predicates()) {
    if (ts_->only_curr(p)) {
      predset_.insert(p);
    }
  }

  // collect variables in bad_
  get_free_symbolic_consts(bad_, vars_in_bad_);
}

// IC3SA specific methods

TermVec IC3SA::symbolic_post_image(size_t i,
                                   UnorderedTermMap & subst,
                                   UnorderedTermMap & last_model_vals)
{
  // TODO use partial_model to handle boolean structure
  //      for now just keeping full boolean structure
  // TODO cache which state updates contain ITEs
  //      and don't even run on the ones which don't

  assert(i >= 1);

  // update last model values
  for (const auto & m : inputvars_at_time_) {
    for (const auto & elem : m) {
      last_model_vals[elem.second] = solver_->get_value(elem.second);
    }
  }

  gen_inputvars_at_time(i - 1);
  const UnorderedTermMap & unconstrained_i_vars = inputvars_at_time_.at(i - 1);
  // update input variables with latest time
  for (const auto & elem : unconstrained_i_vars) {
    subst[elem.first] = elem.second;
    // get the value from the input variable before replacement
    Term val = solver_->get_value(elem.first);
    last_model_vals[elem.second] = val;
  }

  const UnorderedTermMap & state_updates = ts_->state_updates();
  // update state variables value
  for (const auto & sv : ts_->statevars()) {
    if (state_updates.find(sv) != state_updates.end()) {
      subst[sv] = solver_->get_value(sv);
    } else {
      subst[sv] = unconstrained_i_vars.at(sv);
    }
  }

  // TODO: can probably combine this with one of the other loops
  TermVec svs;
  TermVec subst_updates;
  for (const auto & sv : ts_->statevars()) {
    if (state_updates.find(sv) != state_updates.end()) {
      svs.push_back(sv);
      subst_updates.push_back(state_updates.at(sv));
    }
  }

  TermVec replaced_ites = remove_ites_under_model(solver_, subst_updates);
  assert(replaced_ites.size() == svs.size());

  TermVec res;
  // getting sv' = syntactic evaluation of update function
  // and then switching back to sv
  for (size_t i = 0; i < svs.size(); ++i) {
    // TODO: use substitute_terms to save time when copying maps back and forth
    res.push_back(solver_->make_term(
        Equal, svs[i], solver_->substitute(replaced_ites[i], subst)));
  }

  return res;
}

void IC3SA::gen_inputvars_at_time(size_t i)
{
  const UnorderedTermMap & state_updates = ts_->state_updates();
  while (inputvars_at_time_.size() <= i) {
    inputvars_at_time_.push_back(UnorderedTermMap());
    UnorderedTermMap & subst = inputvars_at_time_.back();
    for (const auto & v : ts_->inputvars()) {
      subst[v] = solver_->make_symbol(v->to_string() + "@" + std::to_string(i),
                                      v->get_sort());
    }

    for (auto const & v : ts_->statevars()) {
      if (state_updates.find(v) == state_updates.end()) {
        subst[v] = solver_->make_symbol(
            v->to_string() + "@" + std::to_string(i), v->get_sort());
      }
    }
  }
}

EquivalenceClasses IC3SA::get_equivalence_classes_from_model(
    const UnorderedTermSet & to_keep) const
{
  // assumes the solver state is sat
  EquivalenceClasses ec;
  for (auto elem : term_abstraction_) {
    const Sort & sort = elem.first;
    const UnorderedTermSet & terms = elem.second;

    // TODO figure out if a DisjointSet is a better data structure
    //      will need to keep track of all terms in each partition though

    ec[sort] = std::unordered_map<smt::Term, smt::UnorderedTermSet>();
    std::unordered_map<smt::Term, smt::UnorderedTermSet> & m = ec.at(sort);
    for (auto t : terms) {
      if (in_projection(t, to_keep)) {
        Term val = solver_->get_value(t);
        m[val].insert(t);
      }
    }
  }
  return ec;
}

void IC3SA::construct_partition(const EquivalenceClasses & ec,
                                TermVec & out_cube) const
{
  // now add to the cube expressing this partition
  for (const auto & sortelem : ec) {
    const Sort & sort = sortelem.first;

    // TODO: play around with heuristics for the representative
    //       to add disequalities over
    //       e.g. we're not adding all possible disequalities,
    //       just choosing a representative from each equivalence
    //       class and adding a disequality to encode the distinctness

    //       currently preferring symbol > generic term > value

    // representatives of the different classes of this sort
    TermVec representatives;
    for (const auto & elem : sortelem.second) {
      const Term & val = elem.first;
      assert(val->is_value());

      const UnorderedTermSet & terms = elem.second;
      UnorderedTermSet::const_iterator end = terms.cend();
      UnorderedTermSet::const_iterator it = terms.cbegin();
      Term last = *it;
      it++;

      Term repr = last;
      bool found_repr = false;
      bool repr_val = repr->is_value();

      while (it != end) {
        const Term & term = *(it++);
        assert(last->get_sort() == term->get_sort());
        out_cube.push_back(solver_->make_term(Equal, last, term));
        last = term;

        // TODO: see if a DisjointSet would make this easier
        // update representative for this class
        if (!found_repr) {
          if (term->is_symbolic_const()) {
            repr = term;
            repr_val = false;
            found_repr = true;
          } else if (!term->is_value() && repr_val) {
            repr = term;
            repr_val = false;
          }
        }
      }
    }

    // add disequalities between each pair of representatives from
    // different equivalent classes
    for (size_t i = 0; i < representatives.size(); ++i) {
      for (size_t j = i + 1; j < representatives.size(); ++j) {
        const Term & ti = representatives.at(i);
        const Term & tj = representatives.at(j);
        // should never get the same representative term from different classes
        assert(ti != tj);
        out_cube.push_back(solver_->make_term(Distinct, ti, tj));
      }
    }
  }
}

bool IC3SA::add_to_term_abstraction(const Term & axiom)
{
  SubTermCollector stc(solver_);
  stc.collect_subterms(axiom);

  bool new_terms = false;

  for (const auto & p : stc.get_predicates()) {
    // TODO : figure out if we need to promote all input vars
    //        for this algorithm to work
    //        not sure it's okay to just drop terms containing inputs
    if (ts_->only_curr(p)) {
      new_terms |= predset_.insert(p).second;
    }
  }

  for (const auto & elem : stc.get_subterms()) {
    for (const auto & term : elem.second) {
      // TODO : figure out if we need to promote all input vars
      //        for this algorithm to work
      //        not sure it's okay to just drop terms containing inputs
      if (ts_->only_curr(term)) {
        new_terms |= term_abstraction_[elem.first].insert(term).second;
      }
    }
  }

  return new_terms;
}

void IC3SA::debug_print_equivalence_classes(EquivalenceClasses ec) const
{
  cout << "======== beginning of equivalence classes debug printing" << endl;
  for (const auto & elem : ec) {
    for (const auto & elem2 : elem.second) {
      cout << elem2.first << " : " << elem.first << endl;
      for (const auto & term : elem2.second) {
        cout << "\t" << term << endl;
      }
    }
  }
  cout << "======== end of equivalence classes debug printing" << endl;
}

}  // namespace pono
