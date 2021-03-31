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

namespace cvc4a = ::CVC4::api;
using CVC4SortSet = std::unordered_set<cvc4a::Sort, cvc4a::SortHashFunction>;
using CVC4TermVec = std::vector<cvc4a::Term>;

const unordered_set<PrimOp> bv_ops(
    { Equal, Concat, Extract, BVNot, BVNeg, BVAnd, BVOr, BVXor, BVNand, BVNor,
      BVXnor, BVComp, BVAdd, BVSub, BVMul, BVUdiv, BVSdiv, BVUrem, BVSrem,
      BVSmod, BVShl, BVAshr, BVLshr, BVUlt, BVUle, BVUgt, BVUge, BVSlt, BVSle,
      BVSgt, BVSge, Zero_Extend, Sign_Extend, Repeat, Rotate_Left } );

const unordered_set<PrimOp> relational_ops(
    { And, Or, Xor, Not, Implies, Equal, Distinct, Lt, Le, Gt, Ge,
      BVUlt, BVUle, BVUgt, BVUge, BVSlt, BVSle, BVSgt, BVSge} );

const unordered_map<PrimOp, cvc4a::Kind> to_cvc4_ops(
    { { And, cvc4a::AND },
      { Or, cvc4a::OR },
      { Xor, cvc4a::XOR },
      { Not, cvc4a::NOT },
      { Implies, cvc4a::IMPLIES },
      { Ite, cvc4a::ITE },
      { Equal, cvc4a::EQUAL },
      { Distinct, cvc4a::DISTINCT },
      /* Uninterpreted Functions */
      { Apply, cvc4a::APPLY_UF },
      /* Arithmetic Theories */
      { Plus, cvc4a::PLUS },
      { Minus, cvc4a::MINUS },
      { Negate, cvc4a::UMINUS },
      { Mult, cvc4a::MULT },
      { Div, cvc4a::DIVISION },
      { Lt, cvc4a::LT },
      { Le, cvc4a::LEQ },
      { Gt, cvc4a::GT },
      { Ge, cvc4a::GEQ },
      { Mod, cvc4a::INTS_MODULUS },
      { Abs, cvc4a::ABS },
      { Pow, cvc4a::POW },
      { To_Real, cvc4a::TO_REAL },
      { To_Int, cvc4a::TO_INTEGER },
      { Is_Int, cvc4a::IS_INTEGER },
      /* Fixed Size BitVector Theory */
      { Concat, cvc4a::BITVECTOR_CONCAT },
      // Indexed Op
      { Extract, cvc4a::BITVECTOR_EXTRACT },
      { BVNot, cvc4a::BITVECTOR_NOT },
      { BVNeg, cvc4a::BITVECTOR_NEG },
      { BVAnd, cvc4a::BITVECTOR_AND },
      { BVOr, cvc4a::BITVECTOR_OR },
      { BVXor, cvc4a::BITVECTOR_XOR },
      { BVNand, cvc4a::BITVECTOR_NAND },
      { BVNor, cvc4a::BITVECTOR_NOR },
      { BVXnor, cvc4a::BITVECTOR_XNOR },
      { BVComp, cvc4a::BITVECTOR_COMP },
      { BVAdd, cvc4a::BITVECTOR_PLUS },
      { BVSub, cvc4a::BITVECTOR_SUB },
      { BVMul, cvc4a::BITVECTOR_MULT },
      { BVUdiv, cvc4a::BITVECTOR_UDIV },
      { BVSdiv, cvc4a::BITVECTOR_SDIV },
      { BVUrem, cvc4a::BITVECTOR_UREM },
      { BVSrem, cvc4a::BITVECTOR_SREM },
      { BVSmod, cvc4a::BITVECTOR_SMOD },
      { BVShl, cvc4a::BITVECTOR_SHL },
      { BVAshr, cvc4a::BITVECTOR_ASHR },
      { BVLshr, cvc4a::BITVECTOR_LSHR },
      { BVUlt, cvc4a::BITVECTOR_ULT },
      { BVUle, cvc4a::BITVECTOR_ULE },
      { BVUgt, cvc4a::BITVECTOR_UGT },
      { BVUge, cvc4a::BITVECTOR_UGE },
      { BVSlt, cvc4a::BITVECTOR_SLT },
      { BVSle, cvc4a::BITVECTOR_SLE },
      { BVSgt, cvc4a::BITVECTOR_SGT },
      { BVSge, cvc4a::BITVECTOR_SGE },
      // Indexed Op
      { Zero_Extend, cvc4a::BITVECTOR_ZERO_EXTEND },
      // Indexed Op
      { Sign_Extend, cvc4a::BITVECTOR_SIGN_EXTEND },
      // Indexed Op
      { Repeat, cvc4a::BITVECTOR_REPEAT },
      // Indexed Op
      { Rotate_Left, cvc4a::BITVECTOR_ROTATE_LEFT },
      // Indexed Op
      { Rotate_Right, cvc4a::BITVECTOR_ROTATE_RIGHT },
      // Conversion
      { BV_To_Nat, cvc4a::BITVECTOR_TO_NAT },
      // Indexed Op
      { Int_To_BV, cvc4a::INT_TO_BITVECTOR },
      { Select, cvc4a::SELECT },
      { Store, cvc4a::STORE },
      { Forall, cvc4a::FORALL },
      { Exists, cvc4a::EXISTS },
      // Datatype
      { Apply_Constructor, cvc4a::APPLY_CONSTRUCTOR },
      { Apply_Tester, cvc4a::APPLY_TESTER },
      { Apply_Selector, cvc4a::APPLY_SELECTOR } });

void get_bv_ops_subset(const UnorderedOpSet &in, UnorderedOpSet &out)
{
  for (const auto & o : in) {
    if (bv_ops.find(o.prim_op) != bv_ops.end()) {
      out.insert(o);
    }
  }
}

void collect_values(const Term term, UnorderedTermSet & out)
{
  auto f = [](const smt::Term & t) { return t->is_value(); };
  get_matching_terms(term, out, f);
}

// helper class for generating grammar for CVC4 SyGuS
cvc4a::Grammar cvc4_make_grammar(cvc4a::Solver & cvc4_solver,
                                 const CVC4TermVec & cvc4_boundvars,
                                 const UnorderedOpSet & ops_set,
                                 const vector<cvc4a::Term> * values,
                                 bool all_consts)
{
  // sorts and their terminal constructors (start constructors)
  cvc4a::Sort boolean = cvc4_solver.getBooleanSort();
  cvc4a::Term start_bool = cvc4_solver.mkVar(boolean, "Start");
  vector<cvc4a::Term> start_bvs;

  CVC4SortSet bv_sorts;
  // collect all required sorts
  for (auto cvc4_boundvar : cvc4_boundvars) {
    cvc4a::Sort s = cvc4_boundvar.getSort();
    if (s.isBitVector()) {
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

  // group interesting values according to sorts
  map<cvc4a::Sort, vector<cvc4a::Term>> values_sort_map;
  if (values) {
    for (const auto & v : *values) {
      cvc4a::Sort s = v.getSort();
      auto it = values_sort_map.find(s);
      if (it != values_sort_map.end()) {
        it->second.push_back(v);
      } else {
        values_sort_map[s] = {v};
      }
    }
  }

  if (!ops_set.empty()) {
    vector<cvc4a::Term> g_bound_vars;
    vector<cvc4a::Term> constructs;
    vector<cvc4a::Term> bool_constructs;
    UnorderedOpSet bv_ops_subset;
    get_bv_ops_subset(ops_set, bv_ops_subset);

    for (auto s : start_bvs) {
      g_bound_vars.clear();
      for (auto bound_var : cvc4_boundvars) {
        if (bound_var.getSort() == s.getSort()) {
          g_bound_vars.push_back(bound_var);
        }
      }

      constructs.clear();
      bool_constructs.clear();
      for (auto o : bv_ops_subset) {
        if (get_arity(o.prim_op).first == 1) {
          cout << "UNARY : " << o << endl;
          if (o.prim_op != Extract &&
              o.prim_op != Zero_Extend) {
            constructs.push_back(cvc4_solver.mkTerm(to_cvc4_ops.at(o.prim_op),
                                                    s));
          }
        } else if (get_arity(o.prim_op).first == 2) {
          if (o.prim_op != BVComp &&
              o.prim_op != Concat) {
            cout << "BINARY : " << o << endl;
            if (relational_ops.find(o.prim_op) != relational_ops.end()) {
              bool_constructs.push_back(cvc4_solver.mkTerm(to_cvc4_ops.at(o.prim_op), s, s));
            } else {
              constructs.push_back(cvc4_solver.mkTerm(to_cvc4_ops.at(o.prim_op), s, s));
            }
          } else {
            // TODO: fixme
            continue;
          }
        } else {
          cout << "Unhandled Op : " << o << endl;
          assert(false);
        }
      }

      if (!all_consts) {
        //cvc4a::Term zero = cvc4_solver.mkBitVector(s.getSort().getBVSize(), 0);
        //cvc4a::Term one = cvc4_solver.mkBitVector(s.getSort().getBVSize(), 1);
        //cvc4a::Term min_signed = cvc4_solver.mkBitVector(s.getSort().getBVSize(), pow(2,s.getSort().getBVSize() - 1));
        //constructs.push_back(zero);
        //constructs.push_back(one);
        //constructs.push_back(min_signed);
        auto it = values_sort_map.find(s.getSort());
        if (it != values_sort_map.end()) {
          for (auto v : it->second) {
            constructs.push_back(v);
          }
        }
      } else {
        g.addAnyConstant(s);
      }

      constructs.insert(constructs.end(), g_bound_vars.begin(), g_bound_vars.end());
      g.addRules(s, constructs);
      if (!bool_constructs.empty()) {
        g.addRules(start_bool, bool_constructs);
      }
    }

    // TODO: handle non-bv ops

  } else {
    for (auto s : start_bvs) {
      cvc4a::Term equals = cvc4_solver.mkTerm(cvc4a::EQUAL, s, s);
      cvc4a::Term bvugt = cvc4_solver.mkTerm(cvc4a::BITVECTOR_UGT, s, s);
      g.addRules(start_bool, { bvugt, equals });
    }

    // include bv operations in the grammar
    for (auto s : start_bvs) {
      cvc4a::Term zero = cvc4_solver.mkBitVector(s.getSort().getBVSize(), 0);
      cvc4a::Term one = cvc4_solver.mkBitVector(s.getSort().getBVSize(), 1);
      cvc4a::Term min_signed = cvc4_solver.mkBitVector(s.getSort().getBVSize(), pow(2,s.getSort().getBVSize() - 1));
      cvc4a::Term bvadd = cvc4_solver.mkTerm(cvc4a::BITVECTOR_PLUS, s, s);
      cvc4a::Term bvmul = cvc4_solver.mkTerm(cvc4a::BITVECTOR_MULT, s, s);
      cvc4a::Term bvand = cvc4_solver.mkTerm(cvc4a::BITVECTOR_AND, s, s);
      cvc4a::Term bvcomp = cvc4_solver.mkTerm(cvc4a::BITVECTOR_COMP, s, s);
      cvc4a::Term bvor = cvc4_solver.mkTerm(cvc4a::BITVECTOR_OR, s, s);
      cvc4a::Term bvnot = cvc4_solver.mkTerm(cvc4a::BITVECTOR_NOT, s);
      cvc4a::Term bvneg = cvc4_solver.mkTerm(cvc4a::BITVECTOR_NEG, s);
      vector<cvc4a::Term> g_bound_vars;
      for (auto bound_var : cvc4_boundvars) {
        if (bound_var.getSort() == s.getSort()) {
          g_bound_vars.push_back(bound_var);
        }
      }
      vector<cvc4a::Term> constructs = {bvadd, bvmul,
                                        bvand, bvor, bvnot, bvneg };
      if (!all_consts) {
        constructs.push_back(zero);
        constructs.push_back(one);
        constructs.push_back(min_signed);
      } else {
        g.addAnyConstant(s);
      }
      constructs.insert(
                        constructs.end(), g_bound_vars.begin(), g_bound_vars.end());
      g.addRules(s, constructs);
    }
  }

  return g;
}

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

IC3IA::IC3IA(const Property & p,
             const TransitionSystem & ts,
             const SmtSolver & s,
             PonoOptions opt)
    : super(p, RelationalTransitionSystem(s), s, opt),
      conc_ts_(ts, to_prover_solver_),
      ia_(conc_ts_, ts_, unroller_),
      // only mathsat interpolator supported
      interpolator_(create_interpolating_solver_for(
          SolverEnum::MSAT_INTERPOLATOR, Engine::IC3IA_ENGINE)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      longest_cex_length_(0),
      abs_unroller_(ts_)
{
  engine_ = Engine::IC3IA_ENGINE;
  approx_pregen_ = true;
}

// pure virtual method implementations

IC3Formula IC3IA::get_model_ic3formula() const
{
  TermVec conjuncts;
  conjuncts.reserve(predlbls_.size());
  Term val;
  for (const auto & p : predlbls_) {
    if ((val = solver_->get_value(p)) == solver_true_) {
      conjuncts.push_back(lbl2pred_.at(p));
    } else {
      conjuncts.push_back(solver_->make_term(Not, lbl2pred_.at(p)));
    }
    assert(val->is_value());
  }

  return ic3formula_conjunction(conjuncts);
}

bool IC3IA::ic3formula_check_valid(const IC3Formula & u) const
{
  // check that children are literals
  Term pred;
  Op op;
  for (const auto &c : u.children) {
    if (c->get_sort() != boolsort_) {
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
    if (predset_.find(pred) == predset_.end()) {
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
  if (initialized_) {
    return;
  }

  super::initialize();

  // add all the predicates from init and property to the abstraction
  // NOTE: abstract is called automatically in IC3Base initialize
  UnorderedTermSet preds;
  get_predicates(solver_, conc_ts_.init(), preds, false, false, true);
  size_t num_init_preds = preds.size();
  get_predicates(solver_, bad_, preds, false, false, true);
  size_t num_prop_preds = preds.size() - num_init_preds;
  for (const auto &p : preds) {
    add_predicate(p);
  }
  logger.log(1, "Number predicates found in init: {}", num_init_preds);
  logger.log(1, "Number predicates found in prop: {}", num_prop_preds);
  logger.log(1, "Total number of initial predicates: {}", preds.size());
  // more predicates will be added during refinement
  // these ones are just initial predicates

  // populate cache for existing terms in solver_
  UnorderedTermMap & cache = to_solver_.get_cache();
  Term ns;
  for (auto const&s : ts_.statevars()) {
    // common variables are next states, unless used for refinement in IC3IA
    // then will refer to current state variables after untiming
    // need to cache both
    cache[to_interpolator_.transfer_term(s)] = s;
    ns = ts_.next(s);
    cache[to_interpolator_.transfer_term(ns)] = ns;
  }

  // need to add uninterpreted functions as well
  // first need to find them all
  // NOTE need to use get_free_symbols NOT get_free_symbolic_consts
  // because the latter ignores uninterpreted functions
  UnorderedTermSet free_symbols;
  get_free_symbols(ts_.init(), free_symbols);
  get_free_symbols(ts_.trans(), free_symbols);
  get_free_symbols(bad_, free_symbols);

  for (auto const&s : free_symbols) {
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
  const UnorderedTermSet &bool_symbols = ia_.do_abstraction();

  // don't add boolean symbols that are never used in the system
  // this is an optimization and a fix for some options
  // if using mathsat with bool_model_generation
  // it will fail to get the value of symbols that don't
  // appear in the query
  // thus we don't include those symbols in our cubes
  UnorderedTermSet used_symbols;
  get_free_symbolic_consts(ts_.init(), used_symbols);
  get_free_symbolic_consts(ts_.trans(), used_symbols);
  get_free_symbolic_consts(bad_, used_symbols);

  // add predicates automatically added by ia_
  // to our predset_
  // needed to prevent adding duplicate predicates later
  for (const auto & sym : bool_symbols) {
    assert(sym->is_symbolic_const());
    if (used_symbols.find(sym) != used_symbols.end()) {
      add_predicate(sym);
    }
  }

  assert(ts_.init());  // should be non-null
  assert(ts_.trans());
}

RefineResult IC3IA::refine()
{
  // counterexample trace should have been populated
  assert(cex_.size());
  if (cex_.size() == 1) {
    // if there are no transitions, then this is a concrete CEX
    return REFINE_NONE;
  }

  UnorderedTermSet preds;
  TermVec fresh_preds;

  // HACK added to experiment with CVC4 SyGuS for finding predicates
  if (options_.ic3ia_cvc4_pred_) {
    // first check if the cex trace is spurious
    // this is a bit hacky for now -- should refactor later so that
    // the unrolling is shared with the interpolator approach in the else branch
    // but don't want to disentagle it from interpolation right now
    assert(solver_context_ == 0);
    push_solver_context();

    // NOTE only use abs_unroller in ic3ia-cvc4-pred mode
    // need it to unroll the ts_ in cvc4_find_preds
    // and can't have two unrollers using same solver
    // TODO: revisit the above assumption
    Term bmc_unrolling = abs_unroller_.at_time(conc_ts_.init(), 0);
    Term t;
    bool is_spurious = false;
    Result r;
    for (size_t i = 0; i < cex_.size(); ++i) {
      t = abs_unroller_.at_time(cex_[i], i);
      if (i + 1 < cex_.size()) {
        t = solver_->make_term(And, abs_unroller_.at_time(conc_ts_.trans(), i),
                               t);
      }
      solver_->assert_formula(t);
      r = solver_->check_sat();
      if (r.is_unsat()) {
        is_spurious = true;
        cout << "Shrinking CEX: " << cex_.size() << " -> " << i + 1 << endl;
        cex_.resize(i+1);
        break;
      }
    }

    if (!is_spurious) {
      bmc_unrolling = solver_->make_term(And, bmc_unrolling,
                                         abs_unroller_.at_time(bad_,
                                                               cex_.size() - 1));
      solver_->assert_formula(bmc_unrolling);
      r = solver_->check_sat();
    }

    pop_solver_context();

    assert(!r.is_unknown());
    if (r.is_sat()) {
      // this is a real counterexample
      return REFINE_NONE;
    }
    cvc4_find_preds(cex_, preds);
    fresh_preds.insert(fresh_preds.end(), preds.begin(), preds.end());
  } else {

    size_t cex_length = cex_.size();

    // use interpolator to get predicates
    // remember -- need to transfer between solvers
    assert(interpolator_);

    TermVec formulae;
    for (size_t i = 0; i < cex_length; ++i) {
      // make sure to_solver_ cache is populated with unrolled symbols
      register_symbol_mappings(i);

      Term t = unroller_.at_time(cex_[i], i);
      if (i + 1 < cex_length) {
        t = solver_->make_term(And, t, unroller_.at_time(conc_ts_.trans(), i));
      }
      formulae.push_back(to_interpolator_.transfer_term(t, BOOL));
    }

    TermVec out_interpolants;
    Result r =
      interpolator_->get_sequence_interpolants(formulae, out_interpolants);

    if (r.is_sat()) {
      // this is a real counterexample, so the property is false
      return RefineResult::REFINE_NONE;
    }

    // record the length of this counterexample
    // important to set it here because it's used in register_symbol_mapping
    // to determine if state variables unrolled to a certain length
    // have already been cached in to_solver_
    longest_cex_length_ = cex_length;

    for (auto const&I : out_interpolants) {
      if (!I) {
        assert(
               r.is_unknown());  // should only have null terms if got unknown result
        continue;
      }

      Term solver_I = unroller_.untime(to_solver_.transfer_term(I, BOOL));
      assert(conc_ts_.only_curr(solver_I));
      logger.log(3, "got interpolant: {}", solver_I);
      get_predicates(solver_, solver_I, preds, false, false, true);
    }

    for (auto const&p : preds) {
      if (predset_.find(p) == predset_.end()) {
        // unseen predicate
        fresh_preds.push_back(p);
      }
    }

    if (options_.random_seed_ > 0) {
      shuffle(fresh_preds.begin(),
              fresh_preds.end(),
              default_random_engine(options_.random_seed_));
    }

    // reduce new predicates
    TermVec red_preds;
    if (ia_.reduce_predicates(cex_, fresh_preds, red_preds)) {
      // reduction successful
      logger.log(2,
                 "reduce predicates successful {}/{}",
                 red_preds.size(),
                 fresh_preds.size());
      if (red_preds.size() < fresh_preds.size()) {
        fresh_preds.clear();
        fresh_preds.insert(fresh_preds.end(), red_preds.begin(), red_preds.end());
      }
    } else {
      // should only fail if removed all predicates
      // this can happen when there are uninterpreted functions
      // the unrolling can force incompatible UF interpretations
      // but IC3 (which doesn't unroll) still needs the predicates
      // in this case, just use all the fresh predicates
      assert(red_preds.size() == 0);
      logger.log(2, "reduce predicates FAILED");
    }
  }

  // add all the new predicates
  for (auto const & p : fresh_preds) {
    bool new_pred = add_predicate(p);
    // expect all predicates to be new (e.g. unseen)
    // they were already filtered above
    assert(new_pred);
  }

  logger.log(1, "{} new predicates added by refinement", fresh_preds.size());

  // able to refine the system to rule out this abstract counterexample
  return RefineResult::REFINE_SUCCESS;
}

void IC3IA::reset_solver()
{
  super::reset_solver();

  for (const auto & elem : lbl2pred_) {
    solver_->assert_formula(solver_->make_term(Equal, elem.first, elem.second));
    Term npred = ts_.next(elem.second);
    Term nlbl = label(npred);
    solver_->assert_formula(solver_->make_term(Equal, nlbl, npred));
  }
}

bool IC3IA::is_global_label(const Term & l) const
{
  // all labels used by IC3IA should be globally assumed
  // the assertion will check that this assumption holds though
  assert(super::is_global_label(l) || all_lbls_.find(l) != all_lbls_.end());
  return true;
}

void IC3IA::reabstract()
{
  // don't add boolean symbols that are never used in the system
  // this is an optimization and a fix for some options
  // if using mathsat with bool_model_generation
  // it will fail to get the value of symbols that don't
  // appear in the query
  // thus we don't include those symbols in our cubes
  UnorderedTermSet used_symbols;
  get_free_symbolic_consts(ts_.init(), used_symbols);
  get_free_symbolic_consts(ts_.trans(), used_symbols);
  get_free_symbolic_consts(bad_, used_symbols);

  UnorderedTermSet preds;
  // reset init and trans -- done with calling ia_.do_abstraction
  // then add all boolean constants as (precise) predicates
  for (const auto & p : ia_.do_abstraction()) {
    assert(p->is_symbolic_const());
    if (used_symbols.find(p) != used_symbols.end()) {
      preds.insert(p);
    }
  }

  // predicates from init and bad
  get_predicates(solver_, ts_.init(), preds, false, false, true);
  get_predicates(solver_, bad_, preds, false, false, true);
  // instead of add previously found predicates, we add all the predicates in frame 1
  get_predicates(solver_, get_frame_term(1), preds, false, false, true);

  super::reset_solver();
  if (failed_to_reset_solver_) {
    throw PonoException("IC3IA::reabstract Cannot reabstract because "
                        "the underlying SMT solver doesn't support "
                        "the reset-solver method");
  }
  predset_.clear();
  predlbls_.clear();

  // add predicates
  for (const auto &p : preds) {
    add_predicate(p);
  }
}

bool IC3IA::add_predicate(const Term & pred)
{
  if (predset_.find(pred) != predset_.end()) {
    // don't allow re-adding the same predicate
    return false;
  }

  assert(ts_.only_curr(pred));
  logger.log(2, "adding predicate {}", pred);
  predset_.insert(pred);
  assert(pred->get_sort() == boolsort_);
  assert(pred->is_symbolic_const() || is_predicate(pred, boolsort_));

  Term lbl = label(pred);
  // set the negated label as well
  // can use in either polarity because we add a bi-implication
  labels_[solver_->make_term(Not, pred)] = solver_->make_term(Not, lbl);

  predlbls_.insert(lbl);
  lbl2pred_[lbl] = pred;

  Term npred = ts_.next(pred);
  Term nlbl = label(npred);
  labels_[solver_->make_term(Not, npred)] = solver_->make_term(Not, nlbl);

  if (!pred->is_symbolic_const()) {
    // only need to assert equalities for labels that are distinct
    assert(lbl != pred);
    solver_->assert_formula(solver_->make_term(Equal, lbl, pred));
    solver_->assert_formula(solver_->make_term(Equal, nlbl, npred));

    // only need to modify transition relation for non constants
    // boolean constants will be precise

    // add predicate to abstraction and get the new constraint
    Term predabs_rel = ia_.predicate_refinement(pred);
    static_cast<RelationalTransitionSystem &>(ts_).constrain_trans(predabs_rel);
    // refine the transition relation incrementally
    // by adding a new constraint
    assert(!solver_context_);  // should be at context 0
    solver_->assert_formula(
        solver_->make_term(Implies, trans_label_, predabs_rel));
  }

  // keep track of the labels and different polarities for debugging assertions
  all_lbls_.insert(lbl);
  all_lbls_.insert(solver_->make_term(Not, lbl));
  all_lbls_.insert(nlbl);
  all_lbls_.insert(solver_->make_term(Not, nlbl));

  return true;
}

void IC3IA::register_symbol_mappings(size_t i)
{
  if (i < longest_cex_length_) {
    // these symbols should have already been handled
  }

  UnorderedTermMap & cache = to_solver_.get_cache();
  Term unrolled_sv;
  for (const auto &sv : ts_.statevars()) {
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
  Term abs_trace = abs_unroller_.at_time(ts_.init(), 0);
  Term t;
  for (size_t i = 0; i < cex.size(); ++i) {
    t = cex[i];
    if (i + 1 < cex.size()) {
      t = solver_->make_term(And, t, ts_.trans());
    }
    abs_trace = solver_->make_term(And, abs_trace, abs_unroller_.at_time(t, i));
  }

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
      nv = ts_.next(sv);
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
  size_t num_preds = 1;
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
  logger.log(1, "Looking for {} predicates with CVC4 SyGuS", num_preds);

  bool res = false;

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
  cvc4_solver.setOption("sygus-abort-size", std::to_string(options_.ic3ia_cvc4_pred_size_));

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

  // ops used in the abs_trace
  UnorderedOpSet abs_trace_ops;
  get_ops(abs_trace, abs_trace_ops);
  // constants in the abs_trace
  UnorderedTermSet abs_trace_values;
  collect_values(abs_trace, abs_trace_values);
  vector<cvc4a::Term> abs_trace_values_cvc4;
  for (auto v : abs_trace_values) {
    cvc4a::Term cvc4_val =
      static_pointer_cast<CVC4Term>(to_cvc4_.transfer_term(v))->get_cvc4_term();
    abs_trace_values_cvc4.push_back(cvc4_val);
  }

  cout << "Number of Values : " << abs_trace_values_cvc4.size() << endl;

  // Grammar construction
  cvc4a::Grammar g = cvc4_make_grammar(cvc4_solver, cvc4_boundvars,
                                       abs_trace_ops, NULL,
                                       options_.ic3ia_cvc4_pred_all_consts_);
  cvc4a::Grammar g_with_values =
    cvc4_make_grammar(cvc4_solver, cvc4_boundvars,
                      abs_trace_ops, &abs_trace_values_cvc4,
                      options_.ic3ia_cvc4_pred_all_consts_);

  vector<cvc4a::Term> pred_vec;
  for (size_t n = 0; n < num_preds; ++n) {
    // Create the predicate to search for. Use the grammar
    string pred_name = "P_" + std::to_string(n);
    // cvc4a::Term pred = 
    //   cvc4_solver.synthFun(pred_name, cvc4_boundvars, cvc4_solver.getBooleanSort(), g);
    cvc4a::Term pred;
    switch (n % 3) {
    case 0:
      pred = cvc4_solver.synthFun(pred_name, cvc4_boundvars,
                                  cvc4_solver.getBooleanSort(), g);
      break;
    case 1:
      pred = cvc4_solver.synthFun(pred_name, cvc4_boundvars,
                                  cvc4_solver.getBooleanSort(), g_with_values);
      break;
    default:
      pred = cvc4_solver.synthFun(pred_name, cvc4_boundvars,
                                  cvc4_solver.getBooleanSort());
      break;
    };
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

      cvc4_formula = cvc4_solver.mkTerm(cvc4a::AND,
                                        cvc4_formula,
                                        cvc4_solver.mkTerm(cvc4a::EQUAL, pred_app_vars, pred_app_abs_vars));
    }
  }

  cvc4a::Term constraint = cvc4_solver.mkTerm(cvc4a::NOT, cvc4_formula);

  // use sygus variables rather than ordinary variables.
  std::map<cvc4a::Term, cvc4a::Term> old_to_new;
  std::vector<cvc4a::Term> originals(cvc4_free_vars.begin(),
                                     cvc4_free_vars.end());
  for (cvc4a::Term old_var : originals) {
    cvc4a::Term new_var =
        cvc4_solver.mkSygusVar(old_var.getSort(), old_var.toString() + "_sy");
    old_to_new[old_var] = new_var;
  }
  std::vector<cvc4a::Term> news;
  for (cvc4a::Term old_var : originals) {
    assert(old_to_new.find(old_var) != old_to_new.end());
    news.push_back(old_to_new[old_var]);
  }

  cvc4a::Term sygus_constraint = constraint.substitute(originals, news);
  cvc4_solver.addSygusConstraint(sygus_constraint);
  try {
    res = cvc4_solver.checkSynth().isUnsat();
  }
  catch (cvc4a::CVC4ApiException & e) {
    logger.log(1, "Caught exception from CVC4: {}", e.what());
    return false;
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
    //out_preds.insert(learned_pred);
    get_predicates(solver_, learned_pred, out_preds, false, false, true);
  }

  return res;
}

}  // namespace pono
