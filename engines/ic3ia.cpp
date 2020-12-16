/*********************                                                  */
/*! \file ic3ia.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief IC3 via Implicit Predicate Abstraction (IC3IA) implementation
**        based on
**
**        IC3 Modulo Theories via Implicit Predicate Abstraction
**            -- Alessandro Cimatti, Alberto Griggio,
**               Sergio Mover, Stefano Tonetta
**
**        and the open source implementation:
**
**        https://es-static.fbk.eu/people/griggio/ic3ia/index.html
**
**  within Pono, we are building on the bit-level IC3 instead of directly
**  on IC3Base, because a lot of the functionality is the same
**  In particular, we don't need to override either of the generalization
**  functions. Instead focusing on abstract/refine.
**
**/

#include "engines/ic3ia.h"

#include <random>

#include "smt-switch/cvc4_solver.h"
#include "smt-switch/cvc4_sort.h"
#include "smt-switch/cvc4_term.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"

using namespace smt;
using namespace std;

namespace pono {

IC3IA::IC3IA(Property & p, SolverEnum se, SolverEnum itp_se)
    : super(p, se),
      conc_ts_(property_.transition_system()),
      abs_ts_(solver_),
      ia_(conc_ts_, abs_ts_, unroller_),
      interpolator_(create_interpolating_solver(itp_se)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      longest_cex_length_(0),
      cvc4_(create_solver(smt::SolverEnum::CVC4)),
      to_cvc4_(cvc4_),
      from_cvc4_(solver_)
{
}

IC3IA::IC3IA(Property & p, const SmtSolver & s, SolverEnum itp_se)
    : super(p, s),
      conc_ts_(property_.transition_system()),
      abs_ts_(solver_),
      ia_(conc_ts_, abs_ts_, unroller_),
      interpolator_(create_interpolating_solver(itp_se)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      longest_cex_length_(0),
      cvc4_(create_solver(smt::SolverEnum::CVC4)),
      to_cvc4_(cvc4_),
      from_cvc4_(solver_)
{
}

IC3IA::IC3IA(const PonoOptions & opt,
             Property & p,
             SolverEnum se,
             SolverEnum itp_se)
    : super(opt, p, se),
      conc_ts_(property_.transition_system()),
      abs_ts_(solver_),
      ia_(conc_ts_, abs_ts_, unroller_),
      interpolator_(create_interpolating_solver(itp_se)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      longest_cex_length_(0),
      cvc4_(create_solver(smt::SolverEnum::CVC4)),
      to_cvc4_(cvc4_),
      from_cvc4_(solver_)
{
}

IC3IA::IC3IA(const PonoOptions & opt,
             Property & p,
             const SmtSolver & s,
             SolverEnum itp_se)
    : super(opt, p, s),
      conc_ts_(property_.transition_system()),
      abs_ts_(solver_),
      ia_(conc_ts_, abs_ts_, unroller_),
      interpolator_(create_interpolating_solver(itp_se)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      longest_cex_length_(0),
      cvc4_(create_solver(smt::SolverEnum::CVC4)),
      to_cvc4_(cvc4_),
      from_cvc4_(solver_)
{
}

// pure virtual method implementations

IC3Formula IC3IA::get_model_ic3_formula(TermVec * out_inputs,
                                        TermVec * out_nexts) const
{
  const TermVec & preds = ia_.predicates();
  TermVec conjuncts;
  conjuncts.reserve(preds.size());
  for (auto p : preds) {
    if (solver_->get_value(p) == solver_true_) {
      conjuncts.push_back(p);
    } else {
      conjuncts.push_back(solver_->make_term(Not, p));
    }

    if (out_nexts) {
      Term next_p = ts_->next(p);
      if (solver_->get_value(next_p) == solver_true_) {
        out_nexts->push_back(next_p);
      } else {
        out_nexts->push_back(solver_->make_term(Not, next_p));
      }
    }
  }

  if (out_inputs) {
    for (auto iv : ts_->inputvars()) {
      out_inputs->push_back(
          solver_->make_term(Equal, iv, solver_->get_value(iv)));
    }
  }

  return ic3formula_conjunction(conjuncts);
}

bool IC3IA::ic3formula_check_valid(const IC3Formula & u) const
{
  Sort boolsort = solver_->make_sort(BOOL);
  // check that children are literals
  Term pred;
  Op op;
  for (auto c : u.children) {
    if (c->get_sort() != boolsort) {
      logger.log(3, "ERROR IC3IA IC3Formula contains non-boolean atom: {}", c);
      return false;
    }

    pred = c;
    op = pred->get_op();
    if (op == Not || op == BVNot) {
      pred = *(c->begin());
      assert(pred);
    }

    // expecting either a boolean variable or a predicate
    if (!pred->is_symbolic_const() && predset_.find(pred) == predset_.end()) {
      logger.log(3, "ERROR IC3IA IC3Formula contains unknown atom: {}", pred);
      return false;
    }
  }

  // got through all checks without failing
  return true;
}

void IC3IA::check_ts() const
{
  // basically a No-Op
  // no restrictions except that interpolants must be supported
  // instead of checking explicitly, just let the interpolator throw an
  // exception better than maintaining in two places
}

void IC3IA::initialize()
{
  super::initialize();

  // add all the predicates from init and property to the abstraction
  // NOTE: abstract is called automatically in IC3Base initialize
  UnorderedTermSet preds;
  get_predicates(solver_, ts_->init(), preds, false);
  get_predicates(solver_, bad_, preds, false);
  for (auto p : preds) {
    add_predicate(p);
  }
  // more predicates will be added during refinement
  // these ones are just initial predicates

  // populate cache for existing terms in solver_
  UnorderedTermMap & cache = to_solver_.get_cache();
  Term ns;
  for (auto s : ts_->statevars()) {
    // common variables are next states, unless used for refinement in IC3IA
    // then will refer to current state variables after untiming
    // need to cache both
    cache[to_interpolator_.transfer_term(s)] = s;
    ns = ts_->next(s);
    cache[to_interpolator_.transfer_term(ns)] = ns;
  }

  // need to add uninterpreted functions as well
  // first need to find them all
  // NOTE need to use get_free_symbols NOT get_free_symbolic_consts
  // because the latter ignores uninterpreted functions
  UnorderedTermSet free_symbols;
  get_free_symbols(ts_->init(), free_symbols);
  get_free_symbols(ts_->trans(), free_symbols);
  get_free_symbols(bad_, free_symbols);

  for (auto s : free_symbols) {
    assert(s->is_symbol());
    if (s->is_symbolic_const()) {
      // ignore constants
      continue;
    }
    cache[to_interpolator_.transfer_term(s)] = s;
  }

  // TODO fix generalize_predecessor for ic3ia
  //      might need to override it
  //      behaves a bit differently with both concrete and abstract next state
  //      vars
  if (options_.ic3_pregen_) {
    logger.log(1,
               "WARNING automatically disabling predecessor generalization -- "
               "not supported in IC3IA yet.");
    options_.ic3_pregen_ = false;
  }
}

void IC3IA::abstract()
{
  // main abstraction already done in constructor of ia_
  // just need to set ts_ to the abstraction
  assert(abs_ts_.init());  // should be non-null
  assert(abs_ts_.trans());
  ts_ = &abs_ts_;
}

RefineResult IC3IA::refine()
{
  // recover the counterexample trace
  assert(check_intersects_initial(cex_pg_.target.term));
  TermVec cex({ cex_pg_.target.term });
  ProofGoal tmp = cex_pg_;
  while (tmp.next) {
    tmp = *(tmp.next);
    cex.push_back(tmp.target.term);
    assert(conc_ts_.only_curr(tmp.target.term));
  }

  if (cex.size() == 1) {
    // if there are no transitions, then this is a concrete CEX
    return REFINE_NONE;
  }

  size_t cex_length = cex.size();

  UnorderedTermSet preds;
  // HACK added to experiment with CVC4 SyGuS for finding predicates
  if (options_.ic3ia_cvc4_pred_) {
    cvc4_find_preds(cex, preds);
  } else {
    // use interpolator to get predicates
    // remember -- need to transfer between solvers
    assert(interpolator_);
    Term t = make_and({ conc_ts_.init(), cex[0] });
    Term A = to_interpolator_.transfer_term(unroller_.at_time(t, 0), BOOL);

    TermVec B;
    B.reserve(cex.size() - 1);
    // add to B in reverse order so we can pop_back later
    for (int i = cex.size() - 1; i >= 0; --i) {
      register_symbol_mappings(
          i);  // make sure to_solver_ cache is populated with unrolled symbols
      t = cex[i];
      if (i + 1 < cex.size()) {
        t = solver_->make_term(And, t, conc_ts_.trans());
      }
      B.push_back(
          to_interpolator_.transfer_term(unroller_.at_time(t, i), BOOL));
    }

    // now get interpolants for each transition starting with the first
    bool all_sat = true;
    TermVec interpolants;
    while (B.size()) {
      // Note: have to pass the solver (defaults to solver_)
      Term fullB = make_and(B, interpolator_);
      Term I;
      Result r(smt::UNKNOWN, "unset result");
      try {
        r = interpolator_->get_interpolant(A, fullB, I);
      }
      catch (SmtException & e) {
        logger.log(3, e.what());
      }

      all_sat &= r.is_sat();
      if (r.is_unsat()) {
        Term untimedI = unroller_.untime(to_solver_.transfer_term(I, BOOL));
        logger.log(3, "got interpolant: {}", untimedI);
        interpolants.push_back(untimedI);
      }
      // move next cex time step to A
      // they were added to B in reverse order
      A = interpolator_->make_term(And, A, B.back());
      B.pop_back();
    }

    // record the length of this counterexample
    // important to set it here because it's used in register_symbol_mapping
    // to determine if state variables unrolled to a certain length
    // have already been cached in to_solver_
    longest_cex_length_ = cex.size();

    if (all_sat) {
      // this is a real counterexample, so the property is false
      return RefineResult::REFINE_NONE;
    } else if (!interpolants.size()) {
      logger.log(1,
                 "Interpolation failures...couldn't find any new predicates");
      return RefineResult::REFINE_FAIL;
    } else {
      for (auto I : interpolants) {
        assert(conc_ts_.only_curr(I));
        get_predicates(solver_, I, preds);
      }

      // reduce new predicates
      TermVec preds_vec(preds.begin(), preds.end());
      if (options_.random_seed_ > 0) {
        shuffle(preds_vec.begin(),
                preds_vec.end(),
                default_random_engine(options_.random_seed_));
      }

      TermVec red_preds;
      if (ia_.reduce_predicates(cex, preds_vec, red_preds)) {
        // reduction successful
        preds.clear();
        preds.insert(red_preds.begin(), red_preds.end());
      }
    }
  }

  // add all the new predicates
  bool found_new_preds = false;
  for (auto p : preds) {
    found_new_preds |= add_predicate(p);
  }

  if (!found_new_preds) {
    logger.log(1, "No new predicates found...");
    return RefineResult::REFINE_FAIL;
  }

  // able to refine the system to rule out this abstract counterexample
  return RefineResult::REFINE_SUCCESS;
}

bool IC3IA::add_predicate(const Term & pred)
{
  if (predset_.find(pred) != predset_.end()) {
    // don't allow re-adding the same predicate
    return false;
  }

  assert(ts_->only_curr(pred));
  logger.log(2, "adding predicate {}", pred);
  predset_.insert(pred);
  // add predicate to abstraction and get the new constraint
  Term predabs_rel = ia_.add_predicate(pred);
  // refine the transition relation incrementally
  // by adding a new constraint
  assert(!solver_context_);  // should be at context 0
  solver_->assert_formula(
      solver_->make_term(Implies, trans_label_, predabs_rel));
  return true;
}

void IC3IA::register_symbol_mappings(size_t i)
{
  if (i < longest_cex_length_) {
    // these symbols should have already been handled
  }

  UnorderedTermMap & cache = to_solver_.get_cache();
  Term unrolled_sv;
  for (auto sv : ts_->statevars()) {
    unrolled_sv = unroller_.at_time(sv, i);
    cache[to_interpolator_.transfer_term(unrolled_sv)] = unrolled_sv;
  }
}

bool IC3IA::cvc4_find_preds(const TermVec & cex, UnorderedTermSet & out_preds)
{
  namespace cvc4a = ::CVC4::api;

  // TODO we might need to instantiate a new solver each time? Not sure
  cvc4_->set_opt("incremental", "false");

  // populate from_cvc4_ cache
  UnorderedTermMap & from_cvc4_cache = from_cvc4_.get_cache();
  for (auto sv : conc_ts_.statevars()) {
    from_cvc4_cache[to_cvc4_.transfer_term(sv)] = sv;
  }

  size_t cex_length = cex.size();
  assert(cex_length);

  // unroll the counterexample
  Term solver_formula = unroller_.at_time(abs_ts_.init(), 0);
  Term t;
  for (size_t i = 0; i < cex_length; ++i) {
    t = cex[i];
    if (i + 1 < cex_length) {
      t = solver_->make_term(And, t, abs_ts_.trans());
    }
    solver_formula =
        solver_->make_term(And, solver_formula, unroller_.at_time(t, i));
  }

  // Transfer to a fresh CVC4 solver and get the underlying CVC4 object
  cvc4a::Term cvc4_formula = static_pointer_cast<CVC4Term>(
                                 to_cvc4_.transfer_term(solver_formula, BOOL))
                                 ->get_cvc4_term();

  // get a vector of state variables in the main solver_
  // NOTE: don't forget to use this vector, we're going to rely on the order
  TermVec statevars;
  for (auto sv : conc_ts_.statevars()) {
    statevars.push_back(sv);
  }

  // Retrieve the underlying fresh cvc4 solver
  const cvc4a::Solver & cvc4_solver =
      static_pointer_cast<CVC4Solver>(cvc4_)->get_cvc4_solver();

  // set necessary options for sygus
  cvc4_solver.setOption("lang", "sygus2");
  cvc4_solver.setOption("incremental", "false");

  // create bound variables to use in the synthesized function
  vector<cvc4a::Term> cvc4_statevars;
  cvc4_statevars.reserve(statevars.size());
  vector<cvc4a::Term> cvc4_boundvars;
  cvc4_boundvars.reserve(statevars.size());
  Term transferred_sv;
  for (auto sv : statevars) {
    transferred_sv = to_cvc4_.transfer_term(sv);
    cvc4a::Term cvc4_sv =
        static_pointer_cast<CVC4Term>(transferred_sv)->get_cvc4_term();
    cvc4a::Term cvc4_bv =
        cvc4_solver.mkVar(cvc4_sv.getSort(), cvc4_sv.toString() + "_var");
    cvc4_statevars.push_back(cvc4_sv);
    cvc4_boundvars.push_back(cvc4_bv);
  }

  // Grammr construction
  cvc4a::Sort boolean = cvc4_solver.getBooleanSort();
  cvc4a::Term start_bool = cvc4_solver.mkVar(boolean, "Start");
  std::unordered_set<cvc4a::Sort, cvc4a::SortHashFunction> bv_sorts;
  std::vector<cvc4a::Term> start_bvs;

  for (auto cvc4_boundvar : cvc4_boundvars) {
    cvc4a::Sort s = cvc4_boundvar.getSort();
    bv_sorts.insert(s);
  }
  for (auto s : bv_sorts) {
    cvc4a::Term start_bv = cvc4_solver.mkVar(s, s.toString() + "_start");
    start_bvs.push_back(start_bv);
  }
  vector<cvc4a::Term> starts;
  starts.push_back(start_bool);
  starts.insert(starts.end(), start_bvs.begin(), start_bvs.end());
  cvc4a::Grammar g = cvc4_solver.mkSygusGrammar(cvc4_boundvars, starts);
  
  for (auto s : start_bvs) {
    cvc4a::Term equals = cvc4_solver.mkTerm(cvc4a::EQUAL, s, s);
    cvc4a::Term bvugt = cvc4_solver.mkTerm(cvc4a::BITVECTOR_UGT, s, s);
    g.addRules(start_bool, {equals, bvugt});
  }


  for (auto s : start_bvs) {
    cvc4a::Term zero = cvc4_solver.mkBitVector(s.getSort().getBVSize(), 0);
    cvc4a::Term one = cvc4_solver.mkBitVector(s.getSort().getBVSize(), 1);
    cvc4a::Term bvadd = cvc4_solver.mkTerm(cvc4a::BITVECTOR_PLUS, s, s);
    cvc4a::Term bvmul = cvc4_solver.mkTerm(cvc4a::BITVECTOR_MULT, s, s);
    cvc4a::Term bvand = cvc4_solver.mkTerm(cvc4a::BITVECTOR_AND, s, s);
    cvc4a::Term bvor = cvc4_solver.mkTerm(cvc4a::BITVECTOR_OR, s, s);
    cvc4a::Term bvnot = cvc4_solver.mkTerm(cvc4a::BITVECTOR_NOT, s);
    cvc4a::Term bvneg = cvc4_solver.mkTerm(cvc4a::BITVECTOR_NEG, s);
    vector<cvc4a::Term> g_bound_vars;
    for (auto bound_var : cvc4_boundvars) {
      if (bound_var.getSort() == s.getSort()) {
        g_bound_vars.push_back(bound_var);
      }
    }
    
    vector<cvc4a::Term> constructs = {zero, one, bvadd, bvmul, bvand, bvor, bvnot, bvneg};
    constructs.insert(constructs.end(), g_bound_vars.begin(), g_bound_vars.end());
    g.addRules(s, constructs);
  }

  // Create the predicate to search for
  // TODO: add grammar constraints -- want it to be a predicate
  //       There are also likely heuristic choices in the grammar that will
  //       perform much better
  cvc4a::Term pred =
      cvc4_solver.synthFun("P", cvc4_boundvars, cvc4_solver.getBooleanSort(), g);
//      cvc4_solver.synthFun("P", cvc4_boundvars, cvc4_solver.getBooleanSort());

  // add the implicit predicate abstraction constraints
  // e.g. P(x^@0) <-> P(x@1) /\ P(x^@1) <-> P(x@2) /\ ...
  vector<cvc4a::Term> cvc4_unrolled_next_vars;
  vector<cvc4a::Term> cvc4_unrolled_abstract_vars;
  for (size_t i = 0; i < cex_length; ++i) {
    cvc4_unrolled_next_vars.clear();
    cvc4_unrolled_abstract_vars.clear();

    // CVC4 takes the function as the first child, so insert the pred function
    // first
    cvc4_unrolled_next_vars.push_back(pred);
    cvc4_unrolled_abstract_vars.push_back(pred);

    Term nv, unrolled_nv;
    cvc4a::Term cvc4_unrolled_nv;
    Term abs_nv, unrolled_abs_nv;
    cvc4a::Term cvc4_unrolled_abs_nv;
    for (auto sv : statevars) {
      nv = abs_ts_.next(sv);
      abs_nv = ia_.abstract(nv);

      unrolled_nv = to_cvc4_.transfer_term(unroller_.at_time(nv, i));
      unrolled_abs_nv = to_cvc4_.transfer_term(unroller_.at_time(abs_nv, i));

      cvc4_unrolled_nv =
          static_pointer_cast<CVC4Term>(unrolled_nv)->get_cvc4_term();
      cvc4_unrolled_abs_nv =
          static_pointer_cast<CVC4Term>(unrolled_abs_nv)->get_cvc4_term();

      cvc4_unrolled_next_vars.push_back(cvc4_unrolled_nv);
      cvc4_unrolled_abstract_vars.push_back(cvc4_unrolled_abs_nv);
    }

    cvc4a::Term pred_app_vars =
        cvc4_solver.mkTerm(cvc4a::APPLY_UF, cvc4_unrolled_next_vars);
    cvc4a::Term pred_app_abs_vars =
        cvc4_solver.mkTerm(cvc4a::APPLY_UF, cvc4_unrolled_abstract_vars);

    cvc4_formula = cvc4_solver.mkTerm(
        cvc4a::AND,
        cvc4_formula,
        cvc4_solver.mkTerm(cvc4a::EQUAL, pred_app_vars, pred_app_abs_vars));
  }

  // use sygus variables
  cvc4a::Term constraint = cvc4_solver.mkTerm(cvc4a::NOT, cvc4_formula);
  



  // TODO make sure this is correct -- not sure this makes sense but it should
  // be something like this need to make sure the quantifiers are correct e.g.
  // want to synthesize a predicate such that formula is unsat
  cvc4_solver.addSygusConstraint(constraint);

  bool res = cvc4_solver.checkSynth().isUnsat();
  std::cout << "panda res: " << res << std::endl;
  // for debugging: 
  cvc4_solver.printSynthSolution(std::cout);
  cvc4a::Term pred_solution = cvc4_solver.getSynthSolution(pred);
  std::cout << "panda pred_solution: " << pred_solution << std::endl;
  std::vector<cvc4a::Term> vec;
  vec.push_back(pred_solution);
  vec.insert(vec.end(), std::next(cvc4_unrolled_abstract_vars.begin()), cvc4_unrolled_abstract_vars.end());
  std::cout << "panda size: " << vec.size() << std::endl;
  cvc4a::Term tt = cvc4_solver.mkTerm(cvc4a::APPLY_UF, vec);
  std::cout << "panda application: " << tt << std::endl;

  Term learned_pred;
  // TODO recover the synthesized function and translate it back to a solver_
  // term (learned_pred) and return it

  // will fail here until the TODO above is addressed
  assert(learned_pred);
  // NOTE in future might look for more than one predicate at a time
  out_preds.insert(learned_pred);

  return res;
}

}  // namespace pono
