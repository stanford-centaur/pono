/*********************                                                  */
/*! \file ic3sa.cpp
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
**        Abstraction.
**            -- Aman Goel, Karem Sakallah
**
**
**  within Pono, we are building on the bit-level IC3 instead of directly
**  on IC3Base, because a lot of the functionality is the same
**  In particular, we don't need to override inductive generalization
**
**/

#include "engines/ic3sa.h"

#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "core/prop.h"
#include "core/refineresult.h"
#include "core/rts.h"
#include "core/ts.h"
#include "options/options.h"
#include "smt-switch/result.h"
#include "smt-switch/smt.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "utils/exceptions.h"
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

// main IC3SA implementation

IC3SA::IC3SA(const SafetyProperty & p,
             const TransitionSystem & ts,
             const smt::SmtSolver & solver,
             PonoOptions opt)
    : super(p, RelationalTransitionSystem(solver), solver, opt),
      conc_ts_(ts, to_prover_solver_),
      f_unroller_(conc_ts_, 0, "_AT"),  // zero means pure-functional unrolling
      boolsort_(solver_->make_sort(BOOL)),
      longest_unroll_(0)
{
  // since we passed a fresh RelationalTransitionSystem as the main TS
  // need to point orig_ts_ to the right place
  orig_ts_ = ts;
  engine_ = Engine::IC3SA_ENGINE;
  approx_pregen_ = true;
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

  EquivalenceClasses ec;
  get_equivalence_classes_from_model(ts_.statevars(), ec);
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

  // justify the constraints
  for (const auto & elem : ts_.constraints()) {
    justify_coi(elem.first, all_coi_symbols);
    if (elem.second && ts_.only_curr(elem.first)) {
      justify_coi(ts_.next(elem.first), all_coi_symbols);
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

  EquivalenceClasses ec;
  get_equivalence_classes_from_model(coi_symbols, ec);
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
    if (sk != BOOL && sk != BV && sk != ARRAY) {
      throw PonoException(
          "IC3SA currently only supports bit-vectors and arrays");
    }
  }
}

void IC3SA::abstract()
{
  // need to be able to add path axioms
  // which might require constrain_trans from a relational system
  // ts_ is a relational view of conc_ts_
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

  if (cex_.size() == 1) {
    // TODO if we don't add all terms from init
    // will need to adjust this
    // refine_value doesn't handle it currently, only functional
    // and that could be turned off
    return REFINE_NONE;
  }

  Term learned_lemma;

  RefineResult r;
  bool run_value_refinement = !options_.ic3sa_func_refine_;
  if (options_.ic3sa_func_refine_) {
    // try functional refinement
    logger.log(2, "IC3SA::refine running functional refinement");
    r = ic3sa_refine_functional(learned_lemma);
    // if learned_lemma was simplified to a value, then run value refinement
    if (r == REFINE_SUCCESS) {
      assert(learned_lemma);
      run_value_refinement = learned_lemma->is_value();
    }
  }

  if (run_value_refinement) {
    // possible that simplification reduces the axiom to true
    // it should definitely not be false
    logger.log(2, "IC3SA::refine running value refinement");
    r = ic3sa_refine_value(learned_lemma);
  }

  if (r == REFINE_NONE) {
    assert(!solver_context_);
    return r;
  }

  assert(r == REFINE_SUCCESS);
  assert(learned_lemma);
  assert(!learned_lemma->is_value());

  assert(!ts_.is_functional());
  assert(!solver_context_);
  solver_->assert_formula(
      solver_->make_term(Implies, trans_label_, learned_lemma));
  if (ts_.only_curr(learned_lemma)) {
    ts_.add_constraint(learned_lemma);
    solver_->assert_formula(
        solver_->make_term(Implies, trans_label_, ts_.next(learned_lemma)));
  } else {
    static_cast<RelationalTransitionSystem &>(ts_).constrain_trans(
        learned_lemma);
  }

  logger.log(3, "IC3SA::refine learned axiom {}", learned_lemma);

  // mine for new terms
  // NOTE: can proceed even if there are no new terms because of learned
  // path axiom (learned_lemma)
  UnorderedTermSet new_terms = add_to_term_abstraction(learned_lemma);

  // add symbols from axiom to the projection set permanently
  for (const auto & nt : new_terms) {
    get_free_symbolic_consts(nt, projection_set_);
  }

  if (options_.ic3sa_interp_) {
    assert(interpolator_);
    // TODO: consider only going up until the refined length
    //       currently no way to get it from the refinement functions

    size_t cex_length = cex_.size();
    TermVec formulae;
    for (size_t i = 0; i < cex_length; ++i) {
      // make sure from_interpolator_ cache is populated with unrolled symbols
      register_symbol_mappings(i);

      Term t = unroller_.at_time(cex_[i], i);
      if (i + 1 < cex_length) {
        t = solver_->make_term(And, t, unroller_.at_time(conc_ts_.trans(), i));
      }
      formulae.push_back(to_interpolator_->transfer_term(t, BOOL));
    }

    TermVec out_interpolants;
    Result interp_res =
        interpolator_->get_sequence_interpolants(formulae, out_interpolants);
    assert(!interp_res.is_sat());

    for (auto const & I : out_interpolants) {
      if (!I) {
        // should only have null terms if got unknown result
        assert(interp_res.is_unknown());
        continue;
      }

      Term solver_I =
          unroller_.untime(from_interpolator_->transfer_term(I, BOOL));
      assert(conc_ts_.only_curr(solver_I));
      // TODO look into whether we should add symbols from these terms
      // probably don't need to add symbols from these terms since this is
      // optional
      add_to_term_abstraction(solver_I);
      logger.log(2, "IC3SA::refine interpolant {}", solver_I);
    }

    longest_unroll_ = cex_length;
  }

  assert(!solver_context_);
  return r;
}

RefineResult IC3SA::ic3sa_refine_functional(Term & learned_lemma)
{
  assert(!solver_context_);
  assert(cex_.size());
  // This function will unroll the counterexample trace functionally one step at
  // a time
  // it will introduce fresh symbols for input variables and will keep
  // track of old model values to plug into inputs if an axiom is learned
  Result r;
  UnorderedTermMap last_model_vals;

  assert(!ts_.inputvars().size());
  UnorderedTermSet inputvars;
  // add implicit input variables (states with no update)
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

  if (r.is_sat()) {
    // this is a concrete counterexample
    pop_solver_context();
    assert(!solver_context_);
    return REFINE_NONE;
  }

  assert(r.is_unsat());  // not expecting unknown
  UnorderedTermSet core;
  solver_->get_unsat_assumptions(core);
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

  // This function evaluates each step of the trace using
  // the current state variable values
  // it only keeps one unrolling of the inputs at a time
  // e.g. old unrolled inputs are replaced at the next
  // step by substituting the actual current state variable
  // value and evaluating one transition

  Result r;
  UnorderedTermMap last_model_vals;

  assert(!ts_.inputvars().size());
  UnorderedTermSet inputvars;
  // add implicit input variables (state variables with no update)
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

  if (r.is_sat()) {
    // this is a concrete counterexample
    pop_solver_context();
    assert(!solver_context_);
    return REFINE_NONE;
  }

  assert(r.is_unsat());  // not expecting unknown
  UnorderedTermSet core;
  solver_->get_unsat_assumptions(core);
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

  if (options_.ic3sa_interp_) {
    interpolator_ =
        create_interpolating_solver_for(options_.smt_interpolator_, engine_);
    to_interpolator_ = std::make_unique<TermTranslator>(interpolator_);
    from_interpolator_ = std::make_unique<TermTranslator>(solver_);
  }

  // IC3SA assumes input variables are modeled as state variables
  // with no update function
  // This seems important because otherwise we need to drop terms
  // containing input variables from IC3Formulas
  // It's definitely possible to avoid this (because frames don't
  // need inputs) but in initial tests that doesn't seem to help.
  // It also makes refinement more tricky. There are issues with
  // functional unrolling when the underlying SMT solver simplifies
  // because it might just simplify to an already learned lemma or
  // even true. This happens much more when we don't promote inputs
  // to be state variables and don't allow them in the path axioms
  assert(!ts_.inputvars().size());

  // set up initial term abstraction by getting all subterms

  if (options_.ic3sa_initial_terms_lvl_ <= 2) {
    // for lower options, just add the leaves
    UnorderedTermSet leaves;
    get_leaves(ts_.init(), leaves);
    get_leaves(ts_.trans(), leaves);
    get_leaves(bad_, leaves);
    for (const auto & leaf : leaves) {
      add_to_term_abstraction(leaf);
    }
  }

  if (options_.ic3sa_initial_terms_lvl_ == 1) {
    // also add predicates from bad
    UnorderedTermSet preds;
    get_predicates(solver_, bad_, preds, false, false, true);
    for (const auto & p : preds) {
      add_to_term_abstraction(p);
    }
  } else if (options_.ic3sa_initial_terms_lvl_ == 2) {
    // add predicates from init, trans, and bad
    UnorderedTermSet preds;
    get_predicates(solver_, ts_.init(), preds, false, false, true);
    get_predicates(solver_, ts_.trans(), preds, false, false, true);
    get_predicates(solver_, bad_, preds, false, false, true);
    for (const auto & p : preds) {
      add_to_term_abstraction(p);
    }
  } else if (options_.ic3sa_initial_terms_lvl_ == 3) {
    // only add terms from init and bad
    add_to_term_abstraction(ts_.init());
    add_to_term_abstraction(bad_);
  } else if (options_.ic3sa_initial_terms_lvl_ == 4) {
    // add all terms and predicates
    add_to_term_abstraction(ts_.init());
    add_to_term_abstraction(ts_.trans());
    add_to_term_abstraction(bad_);
  } else {
    // only supports 0-4
    assert(!options_.ic3sa_initial_terms_lvl_);
  }

  Sort boolsort = solver_->make_sort(BOOL);
  // not expecting boolean sorts in term abstraction
  // except for boolector which doesn't distinguish between
  // bit-vectors of size one and booleans
  assert(solver_->get_solver_enum() == BTOR
         || term_abstraction_.find(boolsort) == term_abstraction_.end());
}

// IC3SA specific methods

void IC3SA::get_equivalence_classes_from_model(const UnorderedTermSet & to_keep,
                                               EquivalenceClasses & ec) const
{
  assert(!ec.size());
  // assumes the solver state is sat
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

    //       current priority is: symbol > generic term > value

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
      const Term & ti = representatives.at(i);
      for (size_t j = i + 1; j < representatives.size(); ++j) {
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

UnorderedTermSet IC3SA::add_to_term_abstraction(const Term & term)
{
  UnorderedTermSet new_terms;

  Sort boolsort = solver_->make_sort(BOOL);
  SubTermCollector stc(solver_);
  stc.collect_subterms(term);

  for (const auto & p : stc.get_predicates()) {
    // NOTE might have duplicated predicates wrt equivalences
    //      but ran into issues when I removed it because predicate
    //      values have to be included and currently we're not adding
    //      all possible equalities / disequalities
    //      these are de-duplicated with a set elsewhere
    if (ts_.only_curr(p)) {
      if (predset_.insert(p).second) {
        new_terms.insert(p);
      }
    }
  }

  for (const auto & elem : stc.get_subterms()) {
    for (const auto & term : elem.second) {
      // TODO : figure out if we need to promote all input vars
      //        for this algorithm to work
      //        not sure it's okay to just drop terms containing inputs
      if (ts_.only_curr(term)) {
        if (term_abstraction_[elem.first].insert(term).second) {
          new_terms.insert(term);
        }
      }
    }
  }

  return new_terms;
}

void IC3SA::justify_coi(Term term, UnorderedTermSet & projection)
{
  // expecting to have a satisfiable context
  // and IC3Base only solves at context levels > 0
  assert(solver_context_);

  TermVec to_visit({ term });
  UnorderedTermSet visited;

  Term c;
  while (!to_visit.empty()) {
    c = to_visit.back();
    to_visit.pop_back();

    if (visited.find(c) != visited.end()) {
      // already visited
      continue;
    }
    visited.insert(c);

    Op op = c->get_op();
    Sort sort = c->get_sort();
    if (op == Ite) {
      TermVec children(c->begin(), c->end());
      assert(children.size() == 3);
      to_visit.push_back(children[0]);
      if (solver_->get_value(children[0]) == solver_true_) {
        to_visit.push_back(children[1]);
      } else {
        to_visit.push_back(children[2]);
      }
    } else if (sort == boolsort_
               && is_controlled(op.prim_op, solver_->get_value(c))) {
      to_visit.push_back(get_controlling(c));
    } else {
      for (const auto & cc : c) {
        to_visit.push_back(cc);
      }

      if (ts_.is_next_var(c)) {
        const UnorderedTermMap & state_updates = ts_.state_updates();
        Term curr_c = ts_.curr(c);
        if (state_updates.find(curr_c) != state_updates.end()) {
          to_visit.push_back(state_updates.at(curr_c));
        }
      }

      get_free_symbolic_consts(c, projection);
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

Term IC3SA::get_controlling(Term t) const
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

  Term controlling_term;
  for (const auto & tt : t) {
    if (solver_->get_value(tt) == controlling_val) {
      controlling_term = tt;
      break;
    }
  }
  assert(controlling_term);
  return controlling_term;
}

void IC3SA::register_symbol_mappings(size_t i)
{
  if (i < longest_unroll_) {
    // these symbols should have already been handled
  }

  assert(interpolator_);
  assert(to_interpolator_);
  assert(from_interpolator_);

  UnorderedTermMap & cache = from_interpolator_->get_cache();
  Term unrolled_sv;
  for (const auto & sv : ts_.statevars()) {
    unrolled_sv = unroller_.at_time(sv, i);
    cache[to_interpolator_->transfer_term(unrolled_sv)] = unrolled_sv;
  }
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

bool IC3SA::compute_witness() { return super::compute_witness(conc_ts_); }

}  // namespace pono
