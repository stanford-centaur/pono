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
#include "core/rts.h"
#include "smt-switch/utils.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"
#include "utils/term_walkers.h"

using namespace smt;
using namespace std;

// TODO add option to use interpolants in inductive generalization

namespace pono {

unordered_set<PrimOp> controllable_ops(
    { And,
      Or,
      Implies,
      // include bit-vector versions for boolector
      // will prune out based on sort if
      // not-applicable e.g. for other solvers
      BVOr,
      BVAnd });

// checks if t is an equality literal
// includes BV operators for boolector
// but checks for boolean sort first
bool is_eq_lit(const Term & t, const Sort & boolsort)
{
  const Sort & sort = t->get_sort();
  if (sort != boolsort) {
    return false;
  }

  Op op = t->get_op();
  if (op.is_null()) {
    return false;
  } else if (op == Not || op == BVNot) {
    op = (*(t->begin()))->get_op();
  }

  return (op == Equal || op == BVComp || op == Distinct);
}

// main IC3SA implementation

IC3SA::IC3SA(const Property & p,
             const TransitionSystem & ts,
             const smt::SmtSolver & solver,
             PonoOptions opt)
    : super(p, RelationalTransitionSystem(solver), solver, opt),
      conc_ts_(ts, to_prover_solver_),
      f_unroller_(conc_ts_, 0),  // zero means pure-functional unrolling
      boolsort_(solver_->make_sort(BOOL))
{
  engine_ = Engine::IC3SA_ENGINE;
}

IC3Formula IC3SA::get_model_ic3formula() const
{
  UnorderedTermSet cube_lits;
  // first populate with predicates
  for (const auto & p : predset_) {
    if (solver_->get_value(p) == solver_true_) {
      cube_lits.insert(p);
    } else {
      cube_lits.insert(solver_->make_term(Not, p));
    }
  }

  // TODO make sure that projecting on state variables here makes sense
  EquivalenceClasses ec = get_equivalence_classes_from_model(ts_.statevars());
  construct_partition(ec, cube_lits);
  IC3Formula cube =
      ic3formula_conjunction(TermVec(cube_lits.begin(), cube_lits.end()));
  assert(ic3formula_check_valid(cube));

  return cube;
}

bool IC3SA::ic3formula_check_valid(const IC3Formula & u) const
{
  Sort boolsort = solver_->make_sort(BOOL);
  // check that children are literals
  for (const auto & c : u.children) {
    if (c->get_sort() != boolsort) {
      return false;
    }
  }

  // HACK trouble to accurately identify predicates in boolector
  // because of simplification. Something that starts as a predicate
  // might not end up one, e.g. (= ((_ extract 0 0) x) #b0)
  // just becomes ((_ extract 0 0) x) which we don't consider a predicate
  // TODO possible fix is to just treat all 1-bit vals as predicates in btor
  // TODO check if logging
  if (solver_->get_solver_enum() == BTOR) {
    return true;
  }

  Op op;
  for (auto c : u.children) {
    op = c->get_op();

    // include bvnot for boolector (if we ever include boolector in this check)
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

void IC3SA::predecessor_generalization(size_t i,
                                       const Term & c,
                                       IC3Formula & pred)
{
  UnorderedTermSet all_coi_symbols = projection_set_;
  justify_coi(ts_.next(c), all_coi_symbols);
  assert(all_coi_symbols.size());

  for (const auto & fv : all_coi_symbols) {
    // need to process any constraints that this variable is involved in
    for (const auto & elem : constraint_vars_) {
      const auto & s = elem.second;
      if (s.find(fv) != s.end()) {
        // this variable occurs in this constraint
        // add the constraint
        justify_coi(elem.first, all_coi_symbols);
      }
    }
  }

  // get rid of next-state variables
  UnorderedTermSet coi_symbols;
  for (const auto & tt : all_coi_symbols) {
    if (!ts_.is_next_var(tt)) {
      coi_symbols.insert(tt);
    }
  }
  assert(coi_symbols.size() <= ts_.statevars().size());
  assert(coi_symbols.size());

  logger.log(
      2,
      "IC3SA::generalize_predecessor projecting on {}/{} state variables",
      coi_symbols.size(),
      ts_.statevars().size());


  UnorderedTermSet cube_lits;
  // first populate with predicates
  for (const auto & p : predset_) {
    if (!in_projection(p, coi_symbols)) {
      continue;
    }

    if (solver_->get_value(p) == solver_true_) {
      cube_lits.insert(p);
    } else {
      cube_lits.insert(solver_->make_term(Not, p));
    }
  }

  EquivalenceClasses ec = get_equivalence_classes_from_model(coi_symbols);
  construct_partition(ec, cube_lits);
  pred = ic3formula_conjunction(TermVec(cube_lits.begin(), cube_lits.end()));
  assert(ic3formula_check_valid(pred));
}

void IC3SA::check_ts() const
{
  if (ts_.inputvars().size()) {
    throw PonoException(
        "IC3SA requires all state variables. Try option --promote-inputvars");
  }

  // TODO: add support for arrays

  if (!conc_ts_.is_functional()) {
    throw PonoException("IC3SA requires a functional transition system.");
  }

  for (const auto & sv : ts_.statevars()) {
    SortKind sk = sv->get_sort()->get_sort_kind();
    if (sk != BOOL && sk != BV)
    {
      throw PonoException("IC3SA currently only supports bit-vectors");
    }
  }
  for (const auto & iv : ts_.inputvars()) {
    SortKind sk = iv->get_sort()->get_sort_kind();
    if (sk != BOOL && sk != BV)
    {
      throw PonoException("IC3SA currently only supports bit-vectors");
    }
  }
}

void IC3SA::abstract()
{
  // need to be able to add path axioms
  assert(!ts_.is_functional());
  assert(ts_.solver() == conc_ts_.solver());

  // give this system the same semantics as the functional system
  for (const auto & sv : conc_ts_.statevars()) {
    ts_.add_statevar(sv, conc_ts_.next(sv));
  }
  assert(conc_ts_.inputvars().size() == 0);  // expecting no inputs

  ts_.set_init(conc_ts_.init());

  for (const auto & elem : conc_ts_.state_updates()) {
    ts_.assign_next(elem.first, elem.second);
  }

  for (const auto & elem : conc_ts_.constraints()) {
    ts_.add_constraint(elem.first, elem.second);
  }
}

RefineResult IC3SA::refine()
{
  assert(!solver_context_);
  logger.log(1, "IC3SA: refining a counterexample of length {}", cex_.size());

  Term learned_lemma;

  // try functional refinement
  RefineResult r = ic3sa_refine_functional(learned_lemma);

  if (r == REFINE_NONE) {
    assert(!solver_context_);
    return r;
  }
  assert(r == REFINE_SUCCESS);

  assert(learned_lemma);
  if (learned_lemma->is_value()) {
    // possible that simplification reduces the axiom to true
    // it should definitely not be false
    assert(learned_lemma == solver_true_);
    logger.log(2, "IC3SA::refine falling back on value refinement");
    r = ic3sa_refine_value(learned_lemma);
    assert(learned_lemma);
    assert(!learned_lemma->is_value());
    assert(r == REFINE_SUCCESS);
  }

  assert(!ts_.is_functional());
  assert(!solver_context_);
  solver_->assert_formula(solver_->make_term(Implies, trans_label_, learned_lemma));
  if (ts_.only_curr(learned_lemma)) {
    static_cast<RelationalTransitionSystem &>(ts_).add_constraint(
        learned_lemma);
    solver_->assert_formula(solver_->make_term(Implies,
                                               trans_label_,
                                               ts_.next(learned_lemma)));
  } else {
    static_cast<RelationalTransitionSystem &>(ts_).constrain_trans(
        learned_lemma);
  }

  logger.log(3, "IC3SA::refine learned axiom {}", learned_lemma);

  // add to the projection set permanently
  UnorderedTermSet free_vars;
  get_free_symbolic_consts(learned_lemma, free_vars);
  for (const auto & fv : free_vars) {
    if (ts_.is_curr_var(fv)) {
      projection_set_.insert(fv);
    }
  }

  // mine for new terms
  // NOTE: can proceed even if there are no new terms because of learned
  // path axiom m
  bool new_terms_added = add_to_term_abstraction(learned_lemma);

  // TODO handle refinement failure case if any

  assert(!solver_context_);
  return r;
}

RefineResult IC3SA::ic3sa_refine_functional(Term & learned_lemma)
{
  assert(!solver_context_);
  // This function will unroll the counterexample trace functionally one step at
  // a time it will introduce fresh symbols for input variables it will keep
  // track of old model values to plug into inputs if an axiom is learned
  Result r;
  UnorderedTermMap last_model_vals;

  UnorderedTermSet inputvars = ts_.inputvars();
  // add implicit input variables
  const UnorderedTermMap & state_updates = ts_.state_updates();
  for (const auto & sv : ts_.statevars()) {
    if (state_updates.find(sv) == state_updates.end()) {
      inputvars.insert(sv);
    }
  }

  push_solver_context();

  UnorderedTermSet used_lbls;
  TermVec lbls, assumps;
  Term lbl, unrolled;

  unrolled = f_unroller_.at_time(ts_.init(), 0);
  conjunctive_assumptions(unrolled, used_lbls, lbls, assumps);

  for (size_t i = 0; i < cex_.size(); ++i) {
    unrolled = f_unroller_.at_time(cex_[i], i);
    conjunctive_assumptions(unrolled, used_lbls, lbls, assumps);

    // add constraints
    for (const auto & elem : ts_.constraints()) {
      unrolled = f_unroller_.at_time(elem.first, i);
      conjunctive_assumptions(unrolled, used_lbls, lbls, assumps);
    }

    r = check_sat_assuming(lbls);
    // TODO keep track of model values
    if (r.is_unsat()) {
      break;
    } else {
      assert(r.is_sat());  // not expecting unknown
      // save model values
      Term iv_j;
      for (const auto & iv : inputvars) {
        for (size_t j = 0; j < i; ++j) {
          iv_j = f_unroller_.at_time(iv, j);
          last_model_vals[iv_j] = solver_->get_value(iv_j);
        }
      }
    }
  }
  assert(lbls.size() == assumps.size());

  // TODO check that this correctly handles counterexample case
  if (r.is_sat()) {
    // this is a concrete counterexample
    pop_solver_context();
    assert(!solver_context_);
    return REFINE_NONE;
  }

  assert(r.is_unsat());  // not expecting unknown
  UnorderedTermSet core;
  solver_->get_unsat_core(core);
  assert(core.size());

  TermVec reduced_constraints;
  reduced_constraints.reserve(core.size());
  for (size_t i = 0; i < lbls.size(); ++i) {
    if (core.find(lbls[i]) != core.end()) {
      reduced_constraints.push_back(assumps[i]);
    }
  }
  assert(reduced_constraints.size() == core.size());
  learned_lemma = smart_not(make_and(reduced_constraints));
  learned_lemma = solver_->substitute(learned_lemma, last_model_vals);
  learned_lemma = f_unroller_.untime(learned_lemma);
  assert(ts_.only_curr(learned_lemma));

  pop_solver_context();
  assert(!solver_context_);
  assert(r.is_unsat());
  return REFINE_SUCCESS;
}

RefineResult IC3SA::ic3sa_refine_value(Term & learned_lemma)
{
  assert(!solver_context_);
  assert(cex_.size());

  Result r;
  UnorderedTermMap last_model_vals;

  UnorderedTermSet inputvars = ts_.inputvars();
  // add implicit input variables
  const UnorderedTermMap & state_updates = ts_.state_updates();
  for (const auto & sv : ts_.statevars()) {
    if (state_updates.find(sv) == state_updates.end()) {
      inputvars.insert(sv);
    }
  }

  push_solver_context();

  // due to simplifications can end up with the same terms
  // for the constraints, avoid duplicating labels by keeping
  // track with a set
  UnorderedTermSet used_lbls;
  TermVec lbls, assumps;

  Term p = ts_.init();

  size_t refined_length = 0;
  for (; refined_length + 1 < cex_.size(); ++refined_length) {
    push_solver_context();
    used_lbls.clear();
    lbls.clear();
    assumps.clear();

    conjunctive_assumptions(p, used_lbls, lbls, assumps);
    conjunctive_assumptions(cex_[refined_length], used_lbls, lbls, assumps);
    conjunctive_assumptions(
        ts_.next(cex_[refined_length + 1]), used_lbls, lbls, assumps);
    // if we start with less terms, then need to include trans in the
    // assumptions
    assert_trans_label();

    r = check_sat_assuming(lbls);
    // TODO keep track of model values
    if (r.is_unsat()) {
      break;
    } else {
      assert(r.is_sat());  // not expecting unknown

      UnorderedTermMap subst;

      // save model values
      Term iv_j;
      Term val;
      for (const auto & iv : inputvars) {
        for (size_t j = 0; j <= refined_length; ++j) {
          iv_j = f_unroller_.at_time(iv, j);
          if (j < refined_length) {
            val = solver_->get_value(iv_j);
          } else {
            // wasn't unrolled yet
            val = solver_->get_value(iv);
          }
          last_model_vals[iv_j] = val;
        }
        // unroll inputs
        subst[iv] = f_unroller_.at_time(iv, refined_length);
      }

      // now get the value of the state variables
      for (const auto & elem : state_updates) {
        assert(subst.find(elem.first) == subst.end());
        subst[elem.first] = solver_->get_value(elem.first);
      }

      pop_solver_context();

      // make p the post-state
      p = solver_->make_term(true);
      Term eq;
      for (const auto & elem : state_updates) {
        // TODO look into substitute terms to save time copying over the subst
        // map
        eq = solver_->make_term(
            Equal, elem.first, solver_->substitute(elem.second, subst));
        p = solver_->make_term(And, p, eq);
      }
    }
  }
  assert(lbls.size() == assumps.size());

  // TODO check that this correctly handles counterexample case
  if (r.is_sat()) {
    // this is a concrete counterexample
    pop_solver_context();
    assert(!solver_context_);
    return REFINE_NONE;
  }

  assert(r.is_unsat());  // not expecting unknown
  UnorderedTermSet core;
  solver_->get_unsat_core(core);
  assert(core.size());

  pop_solver_context();
  pop_solver_context();
  assert(!solver_context_);

  TermVec reduced_constraints;
  reduced_constraints.reserve(core.size());
  for (size_t i = 0; i < lbls.size(); ++i) {
    if (core.find(lbls[i]) != core.end()) {
      reduced_constraints.push_back(assumps[i]);
    }
  }
  assert(reduced_constraints.size() == core.size());

  learned_lemma = smart_not(make_and(reduced_constraints));
  learned_lemma = solver_->substitute(learned_lemma, last_model_vals);

  assert(r.is_unsat());
  return REFINE_SUCCESS;
}

void IC3SA::conjunctive_assumptions(const Term & term,
                                    UnorderedTermSet & used_lbls,
                                    TermVec & lbls,
                                    TermVec & assumps)
{
  assert(solver_context_);  // should only add assumptions at non-zero context
  tmp_.clear();
  conjunctive_partition(term, tmp_, true);
  Term lbl;
  for (const auto & tt : tmp_) {
    assert(tt->get_sort() == boolsort_);
    lbl = label(tt);
    if (used_lbls.find(lbl) == used_lbls.end()) {
      used_lbls.insert(lbl);
      lbls.push_back(lbl);
      assumps.push_back(tt);
      solver_->assert_formula(solver_->make_term(Implies, lbl, tt));
    }
  }
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
  // UnorderedTermSet inputvars = ts_.inputvars();
  // for (const auto & iv : inputvars) {
  //   ts_.promote_inputvar(iv);
  // }
  // TODO with change above, don't need to check only_curr everywhere
  // technically should remove those checks (or at least guard with an assert)
  // assert(!ts_.inputvars().size());

  // set up initial term abstraction by getting all subterms
  // TODO consider starting with only a subset -- e.g. variables
  // TODO consider keeping a cache from terms to their free variables
  //      for use in COI

  // TODO make it an option to add ts_.trans()
  // TODO make sure projecting on state variables is right
  // I think we'll always project models onto at least state variables
  // so, we should prune those terms now
  // otherwise we'll do unnecessary iteration over them every time we get a
  // model
  add_to_term_abstraction(ts_.init());
  add_to_term_abstraction(ts_.trans());
  add_to_term_abstraction(bad_);

  Sort boolsort = solver_->make_sort(BOOL);
  // not expecting boolean sorts in term abstraction
  // except for boolector which doesn't distinguish between
  // bit-vectors of size one and booleans
  assert(solver_->get_solver_enum() == BTOR
         || term_abstraction_.find(boolsort) == term_abstraction_.end());

  // collect variables in bad_
  get_free_symbolic_consts(bad_, vars_in_bad_);

  // populate the map used in justify_coi to take constraints into account
  UnorderedTermSet tmp_vars;
  Term next_constraint;
  for (const auto & elem : ts_.constraints()) {
    assert(ts_.no_next(elem.first));
    tmp_vars.clear();
    get_free_symbolic_consts(elem.first, tmp_vars);
    constraint_vars_[elem.first] = tmp_vars;

    if (elem.second) {
      // need to add next state version also
      next_constraint = ts_.next(elem.first);
      tmp_vars.clear();
      get_free_symbolic_consts(next_constraint, tmp_vars);
      constraint_vars_[next_constraint] = tmp_vars;
    }
  }
}

// IC3SA specific methods

EquivalenceClasses IC3SA::get_equivalence_classes_from_model(
    const UnorderedTermSet & to_keep) const
{
  // assumes the solver state is sat
  EquivalenceClasses ec;
  for (const auto & elem : term_abstraction_) {
    const Sort & sort = elem.first;
    const UnorderedTermSet & terms = elem.second;

    // TODO figure out if a DisjointSet is a better data structure
    //      will need to keep track of all terms in each partition though

    ec[sort] = std::unordered_map<smt::Term, smt::UnorderedTermSet>();
    std::unordered_map<smt::Term, smt::UnorderedTermSet> & m = ec.at(sort);
    for (const auto & t : terms) {
      if (in_projection(t, to_keep)) {
        Term val = solver_->get_value(t);
        m[val].insert(t);
      }
    }
  }
  return ec;
}

void IC3SA::construct_partition(const EquivalenceClasses & ec,
                                UnorderedTermSet & out_cube) const
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
    Term lit;
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
        lit = solver_->make_term(Equal, last, term);
        if (!lit->is_value()) {
          // only add if not trivially true
          out_cube.insert(lit);
          assert(solver_->get_value(lit) == solver_true_);
        }
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

      // save the representative for this partition
      representatives.push_back(repr);
    }

    assert(representatives.size() == sortelem.second.size());

    // add disequalities between each pair of representatives from
    // different equivalent classes
    for (size_t i = 0; i < representatives.size(); ++i) {
      for (size_t j = i + 1; j < representatives.size(); ++j) {
        const Term & ti = representatives.at(i);
        const Term & tj = representatives.at(j);
        // should never get the same representative term from different classes
        assert(ti != tj);
        if (ti->is_value() && tj->is_value()) {
          // if there are two values, don't need to express disequality
          continue;
        }
        lit = solver_->make_term(Distinct, ti, tj);
        if (!lit->is_value()) {
          // only add if not trivially true
          out_cube.insert(lit);
          assert(solver_->get_value(lit) == solver_true_);
        }
      }
    }
  }
}

bool IC3SA::add_to_term_abstraction(const Term & term)
{
  Sort boolsort = solver_->make_sort(BOOL);
  SubTermCollector stc(solver_);
  stc.collect_subterms(term);

  bool new_terms = false;

  for (const auto & p : stc.get_predicates()) {
    // TODO : figure out if we need to promote all input vars
    //        for this algorithm to work
    //        not sure it's okay to just drop terms containing inputs
    // NOTE might have duplicated predicates wrt equivalences
    //      but ran into issues when I removed it because predicate
    //      values have to be included and currently we're not adding
    //      all possible equalities / disequalities
    if (ts_.only_curr(p)) {
      new_terms |= predset_.insert(p).second;
    }
  }

  for (const auto & elem : stc.get_subterms()) {
    for (const auto & term : elem.second) {
      // TODO : figure out if we need to promote all input vars
      //        for this algorithm to work
      //        not sure it's okay to just drop terms containing inputs
      if (ts_.only_curr(term)) {
        new_terms |= term_abstraction_[elem.first].insert(term).second;
      }
    }
  }

  return new_terms;
}

void IC3SA::justify_coi(Term c, UnorderedTermSet & projection)
{
  // expecting to have a satisfiable context
  // and IC3Base only solves at context levels > 0
  assert(solver_context_);

  Op op = c->get_op();
  Sort sort = c->get_sort();
  if (op == Ite) {
    TermVec children(c->begin(), c->end());
    assert(children.size() == 3);
    justify_coi(children[0], projection);
    if (solver_->get_value(children[0]) == solver_true_) {
      justify_coi(children[1], projection);
    } else {
      justify_coi(children[2], projection);
    }
  } else if (sort == boolsort_
             && is_controlled(op.prim_op, solver_->get_value(c))) {
    for (const auto & cc : get_controlling(c)) {
      justify_coi(cc, projection);
    }
  } else {
    for (const auto & cc : c) {
      justify_coi(cc, projection);
    }

    if (ts_.is_next_var(c)) {
      const UnorderedTermMap & state_updates = ts_.state_updates();
      Term curr_c = ts_.curr(c);
      if (state_updates.find(curr_c) != state_updates.end()) {
        justify_coi(state_updates.at(curr_c), projection);
      }
    }

    UnorderedTermSet free_vars;
    get_free_symbolic_consts(c, free_vars);
    for (const auto & fv : free_vars) {
      projection.insert(fv);
    }
  }
}

bool IC3SA::is_controlled(PrimOp po, const Term & val) const
{
  assert(val->is_value());

  if (controllable_ops.find(po) == controllable_ops.end()) {
    return false;
  } else if (po == And || po == BVAnd) {
    return val != solver_true_;
  } else if (po == Or || po == BVOr || po == Implies) {
    return val == solver_true_;
  } else {
    throw PonoException("Unhandled case in IC3SA::is_controlled");
  }
}

TermVec IC3SA::get_controlling(Term t) const
{
  assert(solver_context_);

  Op op = t->get_op();
  assert(is_controlled(op.prim_op, solver_->get_value(t)));
  assert(!op.is_null());

  if (op == Implies) {
    // temporarily rewrite t as an Or
    TermVec children(t->begin(), t->end());
    assert(children.size() == 2);
    t = solver_->make_term(Or, smart_not(children[0]), children[1]);
    op = t->get_op();
  }

  assert(op == And || op == BVAnd || op == Or || op == BVOr);

  Term controlling_val = solver_true_;
  if (op == And || op == BVAnd) {
    // controlling val for and is false
    controlling_val = solver_->make_term(false);
  }

  TermVec controlling_terms;
  for (const auto & tt : t) {
    if (solver_->get_value(tt) == controlling_val) {
      controlling_terms.push_back(tt);
    }
  }
  assert(controlling_terms.size());
  return controlling_terms;
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

void IC3SA::debug_print_syntax_abstraction() const
{
  cout << "======== beginning of predicate set debug printing" << endl;
  for (const auto & p : predset_) {
    cout << "\t" << p << endl;
  }
  cout << "======== end of predicate set debug printing" << endl;

  cout << "======== beginning of term abstraction debug printing" << endl;
  for (const auto & sortelem : term_abstraction_) {
    cout << sortelem.first << ":" << endl;
    for (const auto & term : sortelem.second) {
      cout << "\t" << term << endl;
    }
  }
  cout << "======== end of term abstraction debug printing" << endl;
}

void IC3SA::debug_print_cex_trace(const std::vector<IC3Formula> & cex)
{
  cout << "======== beginning of cex trace debug printing" << endl;
  for (size_t i = 0; i < cex.size(); ++i) {
    cout << "c" << i << ": ";
    for (const auto & c : cex[i].children) {
      cout << c << " /\\ ";
    }
    cout << endl;
  }
  cout << "======== end of cex trace debug printing" << endl;
}

}  // namespace pono
