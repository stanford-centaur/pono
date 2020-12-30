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
#include "utils/term_analysis.h"
#include "utils/term_walkers.h"

using namespace smt;
using namespace std;

namespace pono {

bool is_eq_lit(const Term & t, const Sort & boolsort)
{
  if (t->get_sort() != boolsort)
  {
    return false;
  }

  if (t->is_symbolic_const())
  {
    // boolean symbol doesn't need an op
    return true;
  }

  Op op = t->get_op();
  assert(!op.is_null());

  // also include BVNot for boolector
  if (op == Not || op == BVNot)
  {
    op = (*(t->begin()))->get_op();
  }

  if (op != Equal && op != Distinct && op != BVComp)
  {
    return false;
  }

  return true;
}

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
    if (!is_eq_lit(c, boolsort) &&
        !is_lit(c, boolsort) &&
        !is_predicate(c, boolsort, true)) {
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
  assert(solver_context_ == 0);

  // recover the counterexample trace
  assert(check_intersects_initial(cex_pg_->target.term));
  TermVec cex({ cex_pg_->target.term });
  const ProofGoal * tmp = cex_pg_;
  while (tmp->next) {
    tmp = tmp->next;
    cex.push_back(tmp->target.term);
    assert(ts_->only_curr(tmp->target.term));
  }

  size_t cex_length = cex.size();
  assert(cex_length);
  // TODO figure out correct action if cex_length == 1
  if (cex_length == 1) {
    // I suspect these will always be concrete, but need to double check
    throw PonoException("Length one CEX not handled yet");
  }

  // set up initial substitution map
  gen_inputvars_at_time(0);
  UnorderedTermMap subst = inputvars_at_time_.at(0);

  Term p = ts_->init();
  Term c = cex[0];
  Term formula;
  Result r;
  bool refined = false;
  for (size_t i = 1; i < cex_length; ++i) {
    formula = solver_->make_term(And, p, c);
    formula = solver_->make_term(And, formula, ts_->trans());
    formula = solver_->substitute(formula, subst);
    formula = solver_->make_term(And, formula, ts_->next(cex[i]));

    assert(solver_context_ == 0);
    push_solver_context();

    solver_->assert_formula(formula);
    r = solver_->check_sat();
    if (r.is_sat()) {
      p = symbolic_post_image(i, subst, p, c);
      c = cex[i];
    } else {
      refined = true;
      throw PonoException("unsat refinement case not implemented yet");
    }

    pop_solver_context();
  }

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
  throw PonoException("IC3SA::refine NYI");
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

Term IC3SA::symbolic_post_image(size_t i,
                                UnorderedTermMap & subst,
                                const Term & p,
                                const Term & c)
{
  // TODO use partial_model to handle boolean structure
  //      for now just keeping full boolean structure
  // TODO cache which state updates contain ITEs
  //      and don't even run on the ones which don't

  assert(i >= 1);
  Term post_image = solver_->make_term(And, p, c);

  gen_inputvars_at_time(i);
  // update input variables with latest time
  for (const auto & elem : inputvars_at_time_.at(i)) {
    subst[elem.first] = elem.second;
  }

  const UnorderedTermMap & state_updates = ts_->state_updates();
  TermVec svs;
  TermVec subst_updates;
  for (const auto & sv : ts_->statevars()) {
    if (state_updates.find(sv) != state_updates.end()) {
      svs.push_back(sv);
      // NOTE: creating a fresh map with only the relevant substitution
      //       for performance purposes -- otherwise smt-switch has
      //       has to unpack every element of the mapping into
      //       the underlying solver objects
      subst_updates.push_back(
          solver_->substitute(subst.at(sv), { { sv, state_updates.at(sv) } }));
    }
  }

  TermVec replaced_ites = remove_ites_under_model(solver_, subst_updates);
  assert(replaced_ites.size() == svs.size());
  for (size_t i = 0; i < svs.size(); ++i) {
    subst[svs[i]] = replaced_ites[i];
  }

  return solver_->substitute(post_image, subst);
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
      Term repr = val;
      bool found_repr = false;
      bool repr_val = true;

      UnorderedTermSet::const_iterator end = terms.cend();
      UnorderedTermSet::const_iterator it = terms.cbegin();
      Term last = *it;
      it++;

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
