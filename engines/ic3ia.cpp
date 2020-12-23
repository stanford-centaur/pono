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

// includes for the CVC4 pred experiment
#include "smt-switch/cvc4_solver.h"
#include "smt-switch/cvc4_sort.h"
#include "smt-switch/cvc4_term.h"
#include "smt-switch/identity_walker.h"
#include "smt-switch/utils.h"

#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"

using namespace smt;
using namespace std;

namespace pono {

// helper class for translating back from CVC4
// which is not always binarized
class Binarizer : public smt::IdentityWalker
{
  typedef smt::IdentityWalker super;

  // probably this class will inherit from SMT-Switch-Walker
 public:
  Binarizer(smt::SmtSolver & solver) : super(solver, false){};
  ~Binarizer(){};

  smt::Term process(smt::Term t)
  {
    Term res = visit(t);
    return res;
  }

 protected:
  smt::WalkerStepResult visit_term(smt::Term & t) override
  {
    if (!preorder_) {
      Op o = t->get_op();
      PrimOp po = o.prim_op;
      Term res = t;

      TermVec cached_children;
      for (auto tt : t) {
        Term cached_tt;
        assert(in_cache(tt));
        query_cache(tt, cached_tt);
        assert(cached_tt);
        cached_children.push_back(cached_tt);
      }
      switch (po) {
        case BVAnd:
        case BVOr:
        case BVXor:
          // case BVAdd:
        case BVMul:
        case Concat:
          res = solver_->make_term(po, cached_children[0], cached_children[1]);
          for (size_t i = 2; i < cached_children.size(); i++) {
            res = solver_->make_term(po, res, cached_children[i]);
          }
          break;
        default:
          if (!t->is_symbol() && !t->is_value()) {
            res = solver_->make_term(o, cached_children);
          }
          break;
      }

      save_in_cache(t, res);
    }

    return Walker_Continue;
  }
};

IC3IA::IC3IA(Property & p, SolverEnum se, SolverEnum itp_se)
    : super(p, se),
      conc_ts_(property_.transition_system()),
      abs_ts_(solver_),
      ia_(conc_ts_, abs_ts_, unroller_),
      interpolator_(create_interpolating_solver(itp_se)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      longest_cex_length_(0),
      abs_unroller_(abs_ts_, solver_)
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
      abs_unroller_(abs_ts_, solver_)
{
}

IC3IA::IC3IA(Property & p, const SmtSolver & s, SmtSolver itp)
    : super(p, s),
      conc_ts_(property_.transition_system()),
      abs_ts_(solver_),
      ia_(conc_ts_, abs_ts_, unroller_),
      interpolator_(itp),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      longest_cex_length_(0),
      abs_unroller_(abs_ts_, solver_)
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
      abs_unroller_(abs_ts_, solver_)
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
      abs_unroller_(abs_ts_, solver_)
{
}

IC3IA::IC3IA(const PonoOptions & opt,
             Property & p,
             const SmtSolver & s,
             SmtSolver itp)
    : super(opt, p, s),
      conc_ts_(property_.transition_system()),
      abs_ts_(solver_),
      ia_(conc_ts_, abs_ts_, unroller_),
      interpolator_(itp),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      longest_cex_length_(0),
      abs_unroller_(abs_ts_, solver_)
{
}

// pure virtual method implementations

IC3Formula IC3IA::get_model_ic3formula(TermVec * out_inputs,
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

  if (options_.static_coi_) {
    throw PonoException(
        "Abstraction-refinement procedure in IC3IA does not yet work with "
        "static cone-of-influence");
  }
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
    // first check if the cex trace is spurious
    // this is a bit hacky for now -- should refactor later so that
    // the unrolling is shared with the interpolator approach in the else branch
    // but don't want to disentagle it from interpolation right now
    assert(solver_context_ == 0);
    push_solver_context();

    // NOTE only use abs_unroller in ic3ia-cvc4-pred mode
    // need it to unroll the abs_ts_ in cvc4_find_preds
    // and can't have two unrollers using same solver
    Term bmc_unrolling = abs_unroller_.at_time(conc_ts_.init(), 0);
    Term t;
    for (size_t i = 0; i < cex_length; ++i) {
      t = abs_unroller_.at_time(cex[i], i);
      if (i + 1 < cex_length) {
        t = abs_unroller_.at_time(conc_ts_.trans(), i);
      }
      bmc_unrolling = solver_->make_term(And, bmc_unrolling, t);
    }
    bmc_unrolling = solver_->make_term(
        And, bmc_unrolling, abs_unroller_.at_time(bad_, cex_length - 1));
    solver_->assert_formula(bmc_unrolling);
    Result r = solver_->check_sat();

    pop_solver_context();

    assert(!r.is_unknown());
    if (r.is_sat()) {
      // this is a real counterexample
      return REFINE_NONE;
    }
    cvc4_find_preds(cex, preds);
  } else {
    // use interpolator to get predicates
    // remember -- need to transfer between solvers
    assert(interpolator_);
    register_symbol_mappings(0);
    Term t = make_and({ cex[0], conc_ts_.trans() });
    Term A = to_interpolator_.transfer_term(unroller_.at_time(t, 0), BOOL);

    TermVec B;
    B.reserve(cex_length - 1);
    // add to B in reverse order so we can pop_back later
    for (int i = cex_length - 1; i >= 1; --i) {
      // make sure to_solver_ cache is populated with unrolled symbols
      register_symbol_mappings(i);
      t = unroller_.at_time(cex[i], i);
      if (i + 1 < cex_length) {
        t = solver_->make_term(And, t, unroller_.at_time(conc_ts_.trans(), i));
      }
      B.push_back(to_interpolator_.transfer_term(t, BOOL));
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

      // new predicates
      TermVec preds_vec;
      for (auto p : preds) {
        if (predset_.find(p) == predset_.end()) {
          // unseen predicate
          preds_vec.push_back(p);
        }
      }

      if (options_.random_seed_ > 0) {
        shuffle(preds_vec.begin(),
                preds_vec.end(),
                default_random_engine(options_.random_seed_));
      }

      // reduce new predicates
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

  // clear the current proof goals
  // the transitions represented by those backwards reachable traces
  // may not be precise wrt the new predicates
  proof_goals_.clear();

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
  size_t cex_length = cex.size();
  assert(cex_length);

  // unroll the counterexample
  // IMPORTANT NOTE use the abs_unroller_ because unroller_
  // is over conc_ts_
  Term abs_trace = abs_unroller_.at_time(abs_ts_.init(), 0);
  Term t;
  for (size_t i = 0; i < cex_length; ++i) {
    t = cex[i];
    if (i + 1 < cex_length) {
      t = solver_->make_term(And, t, abs_ts_.trans());
    }
    abs_trace = solver_->make_term(And, abs_trace, abs_unroller_.at_time(t, i));
  }
  // redundant but might as well add bad cube
  abs_trace = solver_->make_term(
      And, abs_trace, abs_unroller_.at_time(bad_, cex_length - 1));

  // create unrolled state variables
  // this vector holds pairs where first is next vars
  // and second is abstract variables
  // unrolled at the given time
  vector<pair<TermVec, TermVec>> var_args;

  // free variables in the abs_trace formula
  UnorderedTermSet free_vars;

  // need to use these state variables in order
  // it returns a set but we want to rely on the order
  const UnorderedTermSet & set_statevars = conc_ts_.statevars();
  TermVec statevars(set_statevars.begin(), set_statevars.end());
  size_t num_statevars = statevars.size();

  for (size_t i = 0; i < cex_length; ++i) {
    assert(var_args.size() == i);
    var_args.push_back(make_pair<TermVec, TermVec>(TermVec(), TermVec()));
    TermVec & unrolled_next_vars = var_args[i].first;
    TermVec & unrolled_abs_vars = var_args[i].second;

    unrolled_next_vars.reserve(num_statevars);
    unrolled_abs_vars.reserve(num_statevars);

    Term nv, unrolled_nv;
    Term abs_nv, unrolled_abs_nv;
    for (auto sv : statevars) {
      nv = abs_ts_.next(sv);
      abs_nv = ia_.abstract(nv);

      unrolled_nv = abs_unroller_.at_time(nv, i);
      unrolled_abs_nv = abs_unroller_.at_time(abs_nv, i);

      assert(abs_nv != unrolled_abs_nv);  // check abs vars unrolled corectly

      unrolled_next_vars.push_back(unrolled_nv);
      unrolled_abs_vars.push_back(unrolled_abs_nv);

      // most of these were already added above
      // by traversing cvc4_formula (at smt-switch level)
      // which uses these same variables
      // but this just ensures all variables are in this set
      free_vars.insert(unrolled_nv);
      free_vars.insert(unrolled_abs_nv);

      // just need to get the zero case
      // otherwise covered by next state variable unrolling
      if (i == 0) {
        free_vars.insert(abs_unroller_.at_time(sv, i));
      }
    }

    for (auto iv : conc_ts_.inputvars()) {
      free_vars.insert(abs_unroller_.at_time(iv, i));
    }
  }

  bool found_preds = false;
  size_t num_preds = 0;
  while (!found_preds) {
    found_preds = cvc4_synthesize_preds(
        abs_trace, statevars, var_args, free_vars, num_preds, out_preds);
    num_preds++;
  }

  return found_preds;
}

bool IC3IA::cvc4_synthesize_preds(
    const smt::Term & abs_trace,
    const TermVec & statevars,
    const vector<pair<TermVec, TermVec>> & unrolled_var_args,
    const UnorderedTermSet & free_vars,
    size_t num_preds,
    smt::UnorderedTermSet & out_preds)
{
  bool res = false;
  namespace cvc4a = ::CVC4::api;

  assert(unrolled_var_args.size());

  // HACK
  // hacked in to evaluate CVC4
  // if done for real, should be sure to do this OR the interpolator, not both
  smt::SmtSolver cvc4_ = create_solver(smt::SolverEnum::CVC4);
  smt::TermTranslator to_cvc4_(cvc4_);
  smt::TermTranslator from_cvc4_(solver_);

  // populate from_cvc4_ cache
  UnorderedTermMap & from_cvc4_cache = from_cvc4_.get_cache();
  for (auto sv : conc_ts_.statevars()) {
    from_cvc4_cache[to_cvc4_.transfer_term(sv)] = sv;
  }

  // get the underlying CVC4 objects
  cvc4a::Term cvc4_formula =
      static_pointer_cast<CVC4Term>(to_cvc4_.transfer_term(abs_trace, BOOL))
          ->get_cvc4_term();

  unordered_set<cvc4a::Term, cvc4a::TermHashFunction> cvc4_free_vars;
  for (auto fv : free_vars) {
    cvc4_free_vars.insert(
        static_pointer_cast<CVC4Term>(to_cvc4_.transfer_term(fv))
            ->get_cvc4_term());
  }

  // Retrieve the underlying fresh cvc4 solver
  cvc4a::Solver & cvc4_solver =
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
    cvc4_statevars.push_back(cvc4_sv);
    cvc4a::Term cvc4_bv =
        cvc4_solver.mkVar(cvc4_sv.getSort(), cvc4_sv.toString() + "_var");
    cvc4_boundvars.push_back(cvc4_bv);
  }

  // TODO move this to helper function
  // Grammar construction
  // sorts and their terminal constructors (start constructors)
  cvc4a::Sort boolean = cvc4_solver.getBooleanSort();
  cvc4a::Term start_bool = cvc4_solver.mkVar(boolean, "Start");
  std::unordered_set<cvc4a::Sort, cvc4a::SortHashFunction> bv_sorts;
  std::vector<cvc4a::Term> start_bvs;

  // collect all required sorts
  for (auto cvc4_boundvar : cvc4_boundvars) {
    cvc4a::Sort s = cvc4_boundvar.getSort();
    if (s.isBitVector())
    {
      bv_sorts.insert(s);
    }
  }

  // for each sort, introduce a new constructor for the grammar
  for (auto s : bv_sorts) {
    cvc4a::Term start_bv = cvc4_solver.mkVar(s, s.toString() + "_start");
    start_bvs.push_back(start_bv);
  }

  // merge the Boolean start and the BV start
  vector<cvc4a::Term> starts;
  starts.push_back(start_bool);
  starts.insert(starts.end(), start_bvs.begin(), start_bvs.end());

  // construct the grammar
  cvc4a::Grammar g = cvc4_solver.mkSygusGrammar(cvc4_boundvars, starts);

  for (auto s : start_bvs) {
    cvc4a::Term equals = cvc4_solver.mkTerm(cvc4a::EQUAL, s, s);
    cvc4a::Term bvugt = cvc4_solver.mkTerm(cvc4a::BITVECTOR_UGT, s, s);
    g.addRules(start_bool, {equals, bvugt});
  }

  // include bv operations in the grammar
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

  vector<cvc4a::Term> pred_vec;
  for (size_t n = 0; n < num_preds; ++n) {
    // Create the predicate to search for. Use the grammar
    string pred_name = "P_" + std::to_string(n);
    cvc4a::Term pred =
      cvc4_solver.synthFun(pred_name, cvc4_boundvars, cvc4_solver.getBooleanSort(), g);
    pred_vec.push_back(pred);

    // add the implicit predicate abstraction constraints
    // e.g. P(x^@0) <-> P(x@1) /\ P(x^@1) <-> P(x@2) /\ ...

    vector<cvc4a::Term> cvc4_next_var_args;
    vector<cvc4a::Term> cvc4_abs_var_args;

    for (const auto & var_pair : unrolled_var_args) {
      const TermVec & unrolled_next_vars = var_pair.first;
      const TermVec & unrolled_abs_vars = var_pair.second;

      // CVC4 takes the function as the first child, so insert the
      // pred function first
      cvc4_next_var_args.clear();
      cvc4_next_var_args.push_back(pred);

      cvc4_abs_var_args.clear();
      cvc4_abs_var_args.push_back(pred);

      assert(unrolled_next_vars.size() == unrolled_abs_vars.size());
      Term nv, abs_nv;
      for (size_t i = 0; i < unrolled_next_vars.size(); ++i) {
        nv = to_cvc4_.transfer_term(unrolled_next_vars[i]);
        abs_nv = to_cvc4_.transfer_term(unrolled_abs_vars[i]);

        cvc4_next_var_args.push_back(
            static_pointer_cast<CVC4Term>(nv)->get_cvc4_term());
        cvc4_abs_var_args.push_back(
            static_pointer_cast<CVC4Term>(abs_nv)->get_cvc4_term());
      }

      cvc4a::Term pred_app_vars =
          cvc4_solver.mkTerm(cvc4a::APPLY_UF, cvc4_next_var_args);
      cvc4a::Term pred_app_abs_vars =
          cvc4_solver.mkTerm(cvc4a::APPLY_UF, cvc4_abs_var_args);

      cvc4_formula = cvc4_solver.mkTerm(
                                        cvc4a::AND,
                                        cvc4_formula,
                                        cvc4_solver.mkTerm(cvc4a::EQUAL, pred_app_vars, pred_app_abs_vars));
    }

    cvc4a::Term constraint = cvc4_solver.mkTerm(cvc4a::NOT, cvc4_formula);

    // use sygus variables rather than ordinary variables.
    std::map<cvc4a::Term, cvc4a::Term> old_to_new;
    std::vector<cvc4a::Term> originals(cvc4_free_vars.begin(), cvc4_free_vars.end());
    for (cvc4a::Term old_var : originals) {
      cvc4a::Term new_var = cvc4_solver.mkSygusVar(old_var.getSort(), old_var.toString() + "_sy");
      old_to_new[old_var] = new_var;
    }
    std::vector<cvc4a::Term> news;
    for (cvc4a::Term old_var : originals) {
      assert(old_to_new.find(old_var) != old_to_new.end());
      news.push_back(old_to_new[old_var]);
    }
    cvc4a::Term sygus_constraint = constraint.substitute(originals, news);

    // TODO make sure this is correct -- not sure this makes sense but it should
    // be something like this need to make sure the quantifiers are correct e.g.
    // want to synthesize a predicate such that formula is unsat

    cvc4_solver.addSygusConstraint(sygus_constraint);

    res = cvc4_solver.checkSynth().isUnsat();
  }

  // for debugging:
  cvc4_solver.printSynthSolution(std::cout);

  for (auto pred : pred_vec) {
    cvc4a::Term pred_solution = cvc4_solver.getSynthSolution(pred);

    // instead of applying predicate and simplifying to get rid of the lambda
    // decided to just do substitution
    // this keeps the structure that sygus came up with instead of rewriting it
    // noticed on one example that arithmetic bv operators were replaced with
    // bitwise operators

    vector<cvc4a::Term> pred_solution_children(pred_solution.begin(),
                                               pred_solution.end());
    assert(pred_solution_children.size()
           == 2);  // expecting a bound_var_list and a term
    cvc4a::Term pred_solution_statevars =
      pred_solution_children[1].substitute(cvc4_boundvars, cvc4_statevars);

    Term cvc4_learned_pred = make_shared<CVC4Term>(pred_solution_statevars);

    // Makai: put this back in if it starts failing on term translation
    //        not all backend solvers support n-ary arguments
    // Binarizer binarizer(cvc4_);
    // cvc4_learned_pred = binarizer.process(cvc4_learned_pred);

    Term learned_pred = from_cvc4_.transfer_term(cvc4_learned_pred, BOOL);
    assert(learned_pred);

    // NOTE in future might look for more than one predicate at a time
    out_preds.insert(learned_pred);
  }

  return res;
}

}  // namespace pono
