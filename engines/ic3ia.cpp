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

const unordered_set<cvc4a::Kind> bv_ops({ cvc4a::EQUAL,
                                          cvc4a::BITVECTOR_CONCAT,
                                          cvc4a::BITVECTOR_EXTRACT,
                                          cvc4a::BITVECTOR_NOT,
                                          cvc4a::BITVECTOR_NEG,
                                          cvc4a::BITVECTOR_AND,
                                          cvc4a::BITVECTOR_OR,
                                          cvc4a::BITVECTOR_XOR,
                                          cvc4a::BITVECTOR_NAND,
                                          cvc4a::BITVECTOR_NOR,
                                          cvc4a::BITVECTOR_XNOR,
                                          cvc4a::BITVECTOR_COMP,
                                          cvc4a::BITVECTOR_PLUS,
                                          cvc4a::BITVECTOR_SUB,
                                          cvc4a::BITVECTOR_MULT,
                                          cvc4a::BITVECTOR_UDIV,
                                          cvc4a::BITVECTOR_SDIV,
                                          cvc4a::BITVECTOR_UREM,
                                          cvc4a::BITVECTOR_SREM,
                                          cvc4a::BITVECTOR_SMOD,
                                          cvc4a::BITVECTOR_SHL,
                                          cvc4a::BITVECTOR_ASHR,
                                          cvc4a::BITVECTOR_LSHR,
                                          cvc4a::BITVECTOR_ULT,
                                          cvc4a::BITVECTOR_ULE,
                                          cvc4a::BITVECTOR_UGT,
                                          cvc4a::BITVECTOR_UGE,
                                          cvc4a::BITVECTOR_SLT,
                                          cvc4a::BITVECTOR_SLE,
                                          cvc4a::BITVECTOR_SGT,
                                          cvc4a::BITVECTOR_SGE,
                                          cvc4a::BITVECTOR_ZERO_EXTEND,
                                          cvc4a::BITVECTOR_SIGN_EXTEND,
                                          cvc4a::BITVECTOR_REPEAT,
                                          cvc4a::BITVECTOR_ROTATE_LEFT,
                                          cvc4a::BITVECTOR_ROTATE_RIGHT });

const unordered_set<cvc4a::Kind> relational_ops({
    cvc4a::EQUAL,
    cvc4a::DISTINCT,
    cvc4a::LT,
    cvc4a::LEQ,
    cvc4a::GT,
    cvc4a::GEQ,
    cvc4a::BITVECTOR_ULT,
    cvc4a::BITVECTOR_ULE,
    cvc4a::BITVECTOR_UGT,
    cvc4a::BITVECTOR_UGE,
    cvc4a::BITVECTOR_SLT,
    cvc4a::BITVECTOR_SLE,
    cvc4a::BITVECTOR_SGT,
    cvc4a::BITVECTOR_SGE,
});

const unordered_set<cvc4a::Kind> multisort_ops({ cvc4a::BITVECTOR_EXTRACT,
                                                 cvc4a::BITVECTOR_CONCAT,
                                                 cvc4a::BITVECTOR_ZERO_EXTEND,
                                                 cvc4a::BITVECTOR_COMP });

const unordered_set<cvc4a::Kind> unary_ops({ cvc4a::BITVECTOR_NEG,
                                             cvc4a::BITVECTOR_NOT,
                                             cvc4a::BITVECTOR_EXTRACT,
                                             cvc4a::BITVECTOR_ZERO_EXTEND,
                                             cvc4a::BITVECTOR_SIGN_EXTEND,
                                             cvc4a::UMINUS });

const unordered_set<cvc4a::Kind> bool_ops(
    { cvc4a::AND, cvc4a::OR, cvc4a::XOR, cvc4a::NOT, cvc4a::IMPLIES, cvc4a::ITE });

// Helpers for CVC4 SyGuS Predicate Search
// should eventually be moved elsewhere

bool cvc4_term_is_value(const cvc4a::Term & term)
{
  cvc4a::Kind k = term.getKind();
  return ((k == cvc4a::CONST_BOOLEAN) || (k == cvc4a::CONST_BITVECTOR)
          || (k == cvc4a::CONST_RATIONAL) || (k == cvc4a::CONST_FLOATINGPOINT)
          || (k == cvc4a::CONST_ROUNDINGMODE) || (k == cvc4a::CONST_STRING)
          || (k == cvc4a::CONST_ARRAY));
}

using CVC4OpSignatures =
    unordered_map<cvc4a::Op,
                  unordered_set<cvc4a::Sort, cvc4a::SortHashFunction>,
                  cvc4a::OpHashFunction>;
using CVC4ValueMap =
    unordered_map<cvc4a::Sort,
                  unordered_set<cvc4a::Term, cvc4a::TermHashFunction>,
                  cvc4a::SortHashFunction>;

/** \class CVC4GrammarSeed
 *  \brief A class for seeding a SyGuS grammar with terms
 *
 *  Used to store ops and values used in an abstract trace
 *  More specifically, it keeps track not only of the operators,
 *  but also which sorts they're applied to
 */
class CVC4GrammarSeed
{
 public:
  CVC4GrammarSeed(cvc4a::Solver & solver) : solver_(solver), num_values_(0) {}

  void scan(cvc4a::Term term)
  {
    CVC4TermVec to_visit({ term });
    unordered_set<cvc4a::Term, cvc4a::TermHashFunction> visited;
    unordered_set<cvc4a::Term, cvc4a::TermHashFunction> contains_free_vars;
    cvc4a::Term t;
    while (!to_visit.empty()) {
      t = to_visit.back();
      cvc4a::Sort sort = t.getSort();
      to_visit.pop_back();
      if (visited.find(t) != visited.end()) {
        for (const auto & tt : t) {
          if (contains_free_vars.find(tt) != contains_free_vars.end()) {
            // mark this term as containing free variables
            contains_free_vars.insert(t);
            break;
          }
        }

        // HACK want to avoid adding integer sort
        // if it's only for values (e.g. 2)
        // but then all variables are reals
        // NOTE: can't just use the CVC4 equivalent of is_value
        // because it would return false for (- 1)
        // since it has an operator
        if (!sort.isInteger()
            || contains_free_vars.find(t) != contains_free_vars.end()) {
          all_sorts_.insert(sort);
        }

        continue;
      } else {
        visited.insert(t);
        to_visit.push_back(t);
        to_visit.insert(to_visit.end(), t.begin(), t.end());

        if (cvc4_term_is_value(t)) {
          value_map_[sort].insert(t);
          num_values_++;
        } else if (t.hasOp()) {
          cvc4a::Op op = t.getOp();
          // easiest way to store signature is as a function sort
          vector<cvc4a::Sort> sort_vec;
          for (const auto & tt : t) {
            sort_vec.push_back(tt.getSort());
          }

          cvc4a::Sort ufsort = solver_.mkFunctionSort(sort_vec, sort);
          assert(!op.isNull());
          op_map_[op].insert(ufsort);
        } else if (t.getKind() == cvc4a::CONSTANT) {
          contains_free_vars.insert(t);
        }
      }
    }
  }

  /** Getter for operator map
   *  @return a map from Ops to signatures (stored as function sorts)
   *  e.g. bvadd might have been applied to (_ BitVec 8) x (_ BitVec 8) -> (_
   * BitVec 8) and (_ BitVec 6) x (_ BitVec 6) -> (_ BitVec 6)
   */
  const CVC4OpSignatures & get_op_map() const { return op_map_; }

  /** Getter for value map
   *  @return map from sorts to a set of values used of that sort
   */
  const CVC4ValueMap & get_value_map() const { return value_map_; }

  const unordered_set<cvc4a::Sort, cvc4a::SortHashFunction> & get_all_sorts()
      const
  {
    return all_sorts_;
  }

  size_t num_values() const { return num_values_; }

 protected:
  cvc4a::Solver & solver_;
  CVC4OpSignatures op_map_;
  CVC4ValueMap value_map_;
  unordered_set<cvc4a::Sort, cvc4a::SortHashFunction> all_sorts_;
  size_t num_values_;
};

// helper class for generating grammar for CVC4 SyGuS
cvc4a::Grammar cvc4_make_grammar(
    cvc4a::Solver & cvc4_solver,
    const CVC4TermVec & cvc4_boundvars,
    const CVC4GrammarSeed * gs,
    unordered_map<cvc4a::Sort,
                  unordered_set<cvc4a::Term, cvc4a::TermHashFunction>,
                  cvc4a::SortHashFunction> & cvc4_max_terms,
    // 0 - no values
    // 1 - values from CVC4GrammarSeed
    // 2 - all values
    size_t values,
    bool all_sorts)
{
  // sorts and their terminal constructors (start constructors)
  cvc4a::Sort boolean = cvc4_solver.getBooleanSort();
  cvc4a::Term start_bool = cvc4_solver.mkVar(boolean, "Start");
  vector<cvc4a::Term> start_terms;

  unordered_map<cvc4a::Sort, cvc4a::Term, cvc4a::SortHashFunction> sort2start(
      { { boolean, start_bool } });
  for (const auto & cvc4_bv : cvc4_boundvars)
  {
    cvc4a::Sort sort = cvc4_bv.getSort();
    // TODO: remove this limitation
    if (!sort.isBitVector() && !sort.isBoolean() && !sort.isReal()
        && !sort.isInteger()) {
      cout << "Skipping unsupported sort " << sort << endl;
      continue;
    }

    if (sort2start.find(sort) == sort2start.end())
    {
      cvc4a::Term start_term = cvc4_solver.mkVar(sort, sort.toString() + "_start");
      sort2start[sort] = start_term;
      start_terms.push_back(start_term);
    }
  }

  unordered_map<cvc4a::Term, vector<cvc4a::Term>, cvc4a::TermHashFunction>
      constructs;

  // collect variables by sort and add to constructs
  unordered_map<cvc4a::Sort, vector<cvc4a::Term>, cvc4a::SortHashFunction>
      sort2boundvars;
  for (const auto & bvar : cvc4_boundvars) {
    cvc4a::Sort sort = bvar.getSort();
    sort2boundvars[sort].push_back(bvar);
    constructs[sort2start.at(sort)].push_back(bvar);
  }

  if (gs) {
    const CVC4OpSignatures & ops_map = gs->get_op_map();
    const CVC4ValueMap & values_map = gs->get_value_map();

    if (all_sorts) {
      // add all sorts to start terms
      for (const auto & sort : gs->get_all_sorts()) {
        // treat booleans specially, only included as result of relational
        // operations
        if (sort2start.find(sort) == sort2start.end() && !sort.isBoolean()) {
          cvc4a::Term start_term =
              cvc4_solver.mkVar(sort, sort.toString() + "_start");
          sort2start[sort] = start_term;
          start_terms.push_back(start_term);
        }
      }
    }

    // separate operators
    unordered_map<SortKind, unordered_set<cvc4a::Kind>> reg_ops;
    unordered_map<SortKind, unordered_set<cvc4a::Kind>> rel_ops;
    unordered_map<SortKind, unordered_set<cvc4a::Op, cvc4a::OpHashFunction>>
        ms_ops;

    assert(ops_map.size());
    for (const auto & opelem : ops_map) {
      cvc4a::Op op = opelem.first;
      cvc4a::Kind po = op.getKind();
      // for now, just bv or arithmetic (using REAL for both real and integer)
      SortKind sk = bv_ops.find(po) == bv_ops.end() ? REAL : BV;
      if (multisort_ops.find(po) != multisort_ops.end()) {
        ms_ops[sk].insert(op);
      } else if (relational_ops.find(po) != relational_ops.end()) {
        assert(!op.isIndexed());
        rel_ops[sk].insert(po);
      } else if (bool_ops.find(po) != bool_ops.end()) {
        // skip boolean operators -- looking for predicates
        assert(!op.isIndexed());
      } else {
        assert(!op.isIndexed());
        reg_ops[sk].insert(po);
      }
    }
    assert(ms_ops.size() + rel_ops.size() + reg_ops.size());

    // regular and relational operators
    for (const auto & s : start_terms) {
      cvc4a::Sort sort = s.getSort();
      SortKind sk;
      if (sort.isBitVector()) {
        sk = BV;
      } else if (sort.isReal() || sort.isInteger()) {
        // using real for both integers and reals
        sk = REAL;
      } else if (sort.isBoolean()) {
        // only looking for predicates
        // nothing to do with a boolean
        // shouldn't be included in start terms
        assert(false);
      } else {
        // TODO handle INT/REAL (note using REAL for both)
        logger.log(0,
                   "WARNING IC3IA CVC4 predicates unhandled sort: {}",
                   s.getSort().toString());
      }

      // regular
      for (const auto & po : reg_ops[sk]) {
        if (unary_ops.find(po) == unary_ops.end()) {
          constructs[s].push_back(cvc4_solver.mkTerm(po, s, s));
        } else {
          constructs[s].push_back(cvc4_solver.mkTerm(po, s));
        }
      }

      // relational
      for (const auto & po : rel_ops[sk]) {
        constructs[start_bool].push_back(cvc4_solver.mkTerm(po, s, s));
      }

      // add values
      if (values == 1) {
        auto it = values_map.find(s.getSort());
        if (it != values_map.end()) {
          for (const auto & val : it->second) {
            constructs[s].push_back(val);
          }
        }
      }
    }

    // handle multi-sort operators
    // currently only for bit-vectors
    for (SortKind sk : { BV }) {
      for (const auto & op : ms_ops[sk]) {
        for (const auto & signature : ops_map.at(op)) {
          assert(signature.isFunction());
          vector<cvc4a::Term> args;
          bool failed = false;
          for (const auto & sort : signature.getFunctionDomainSorts()) {
            // TODO might not have start terms for all sorts
            //      need to create a start term for each sort in the problem
            if (sort2start.find(sort) == sort2start.end()) {
              // no start term for this sort
              failed = true;
              logger.log(
                  3,
                  "IC3IA SyGuS pred: skipping op {} because missing sort {}",
                  op.toString(),
                  sort.toString());
              break;
            }
            args.push_back(sort2start.at(sort));
          }
          if (failed) {
            break;
          }

          cvc4a::Sort codomain_sort = signature.getFunctionCodomainSort();
          if (sort2start.find(codomain_sort) != sort2start.end()) {
            cvc4a::Term return_start_term = sort2start.at(codomain_sort);
            constructs[return_start_term].push_back(
                cvc4_solver.mkTerm(op, args));
          } else {
            logger.log(3,
                       "Skipping op {} because missing return sort {}",
                       op.toString(),
                       codomain_sort.toString());
          }
        }
      }
    }

  } else {
    for (auto s : start_terms) {
      cvc4a::Term equals = cvc4_solver.mkTerm(cvc4a::EQUAL, s, s);
      cvc4a::Term bvugt = cvc4_solver.mkTerm(cvc4a::BITVECTOR_UGT, s, s);
      constructs[start_bool] = { bvugt, equals };
    }

    // include bv operations in the grammar
    for (auto s : start_terms) {
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
      constructs[s] = { bvadd, bvmul, bvand, bvor, bvnot, bvneg };
      if (values == 0) {
        // no constants in this case
        ;
      } else if (values == 1) {
        constructs[s].push_back(zero);
        constructs[s].push_back(one);
        constructs[s].push_back(min_signed);
      } else {
        assert(values == 2);
      }
    }

    // TODO: non-bv ops
  }

  // add max terms
  for (const auto & start_term : start_terms) {
    cvc4a::Sort sort = start_term.getSort();
    for (const auto & t : cvc4_max_terms[sort]) {
      constructs[start_term].push_back(t);
    }
  }

  for (const auto & start_term : start_terms) {
    cvc4a::Sort sort = start_term.getSort();
    if (sort2boundvars.find(sort) == sort2boundvars.end()) {
      // there should be some kind of start (e.g. value or constant)
      // for every start term
      // it's possible there will be no constants if the only terms
      // of a given sort contain input variables
      if (!sort.isBitVector()) {
        throw PonoException("Unhandled empty start term of sort: "
                            + sort.toString());
      }

      cvc4a::Term zero = cvc4_solver.mkBitVector(sort.getBVSize(), 0);
      cvc4a::Term one = cvc4_solver.mkBitVector(sort.getBVSize(), 1);
      constructs[start_term].push_back(zero);
      constructs[start_term].push_back(one);
    }
  }

  // construct the grammar
  vector<cvc4a::Term> starts({ start_bool });
  starts.reserve(start_terms.size() + 1);
  starts.insert(starts.end(), start_terms.begin(), start_terms.end());
  cvc4a::Grammar g = cvc4_solver.mkSygusGrammar(cvc4_boundvars, starts);

  for (const auto & elem : constructs) {
    const cvc4a::Term & start_term = elem.first;
    const vector<cvc4a::Term> & rules = elem.second;
    assert(rules.size());

    // add rules for this nonterminal start term
    g.addRules(start_term, rules);

    if (values == 2) {
      g.addAnyConstant(start_term);
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
      // TODO: when cleaning this up, should have only one unroller
      ia_(conc_ts_, ts_, opt.ic3ia_cvc4_pred_ ? abs_unroller_ : unroller_),
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

  // collect all sort kinds
  auto varsets = { ts_.statevars(), ts_.inputvars() };
  for (const auto & vs : varsets) {
    for (const auto & v : vs) {
      all_sort_kinds_.insert(v->get_sort()->get_sort_kind());
    }
  }

  // find candidate predicates from trans
  UnorderedTermSet pred_cand_set;
  if (conc_ts_.is_functional()) {
    for (const auto & elem : conc_ts_.state_updates()) {
      get_predicates(solver_, elem.second, pred_cand_set, false, true, true);
    }
    for (const auto & elem : conc_ts_.constraints()) {
      get_predicates(solver_, elem.first, pred_cand_set, false, true, true);
    }
  } else {
    get_predicates(solver_, conc_ts_.trans(), pred_cand_set, false, true, true);
  }

  for (const auto & p : predset_) {
    pred_cand_set.erase(p);
  }

  assert(!pred_candidates_.size());
  for (const auto & p : pred_cand_set) {
    if (conc_ts_.only_curr(p)) {
      pred_candidates_.push_back(p);
    }
  }
  logger.log(1,
             "Found {} predicate candidate(s) in transition relation",
             pred_candidates_.size());

  // gather max-terms (for ic3ia-cvc4-pred)
  if (options_.ic3ia_cvc4_pred_maxterms_) {
    UnorderedTermSet all_preds = preds;
    SolverEnum solver_enum = solver_->get_solver_enum();
    if (conc_ts_.is_functional()) {
      for (const auto & elem : conc_ts_.state_updates()) {
        Sort esort = elem.second->get_sort();
        // constants and values get added anyway, don't count as max term
        if (esort != boolsort_ && conc_ts_.only_curr(elem.second)
            && !elem.second->is_value() && !elem.second->is_symbol()) {
          max_terms_.insert(elem.second);
        }
      }
      for (const auto & elem : conc_ts_.constraints()) {
        get_predicates(solver_, elem.first, all_preds, false, false, false);
      }

      for (Term p : all_preds) {
        if (p->is_symbolic_const()) {
          continue;
        }

        Sort tt_sort;
        for (Term tt : p) {
          tt_sort = tt->get_sort();
          assert(solver_enum == BTOR || tt_sort != boolsort_);
          // constants and values get added anyway, don't count as max term
          if (conc_ts_.only_curr(tt) && !tt->is_value() && !tt->is_symbol()) {
            // TODO: consider calling promote-inputvars always to avoid this
            // issue
            max_terms_.insert(tt);
          }
        }
      }
    } else {
      get_predicates(solver_, conc_ts_.trans(), all_preds, false, false, false);

      for (Term p : all_preds) {
        if (p->is_symbolic_const()) {
          continue;
        }

        // recurse until hit largest only_curr term
        TermVec to_process(p->begin(), p->end());
        while (to_process.size()) {
          Term tt = to_process.back();
          to_process.pop_back();
          // constants and values get added anyway, don't count as max term
          if (conc_ts_.only_curr(tt) && !tt->is_value() && !tt->is_symbol()) {
            Sort s = tt->get_sort();
            if (s == boolsort_) {
              continue;
            }
            max_terms_.insert(tt);
          } else {
            for (Term c : tt) {
              to_process.push_back(c);
            }
          }
        }
      }
    }
  }

  logger.log(1, "IC3IA: Got {} max terms", max_terms_.size());
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
        logger.log(1, "Shrinking CEX: {} -> {}", cex_.size(), i + 1);
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

    // try identifying a predicate from known predicate candidates first
    // (found in the trans of the concrete system)
    TermVec pred_vec;
    if (pred_candidates_.size()
        && ia_.reduce_predicates(cex_, pred_candidates_, pred_vec)) {
      logger.log(
          1, "Found {} sufficient predicates from candidates", pred_vec.size());
      preds.insert(pred_vec.begin(), pred_vec.end());
    } else {
      assert(!preds.size());
      // otherwise fall back on SyGuS
      cvc4_find_preds(cex_, preds);
    }

    for (auto const & p : preds) {
      if (predset_.find(p) == predset_.end()) {
        // unseen predicate
        fresh_preds.push_back(p);
      }
    }
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
  int pred_size = options_.ic3ia_cvc4_pred_size_;
  pred_size += (num_preds-1)/3; // increase the size periodically
  cvc4_solver.setOption("sygus-abort-size", std::to_string(pred_size));
  cvc4_solver.setOption("sygus-active-gen", "enum");

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

  unordered_map<cvc4a::Sort,
                unordered_set<cvc4a::Term, cvc4a::TermHashFunction>,
                cvc4a::SortHashFunction>
      cvc4_max_terms;
  for (const auto & t : max_terms_) {
    cvc4a::Term cvc4_t =
        static_pointer_cast<CVC4Term>(to_cvc4_.transfer_term(t))
            ->get_cvc4_term();
    // need to substitute with bound vars
    cvc4_max_terms[cvc4_t.getSort()].insert(
        cvc4_t.substitute(cvc4_statevars, cvc4_boundvars));
  }

  CVC4GrammarSeed gs(cvc4_solver);
  Term transferred_trace = to_cvc4_.transfer_term(abs_trace, BOOL);
  gs.scan(static_pointer_cast<CVC4Term>(transferred_trace)->get_cvc4_term());

  logger.log(1, "Number of Values : {}", gs.num_values());

  // Grammar construction
  cvc4a::Grammar g =
      cvc4_make_grammar(cvc4_solver,
                        cvc4_boundvars,
                        &gs,
                        cvc4_max_terms,
                        options_.ic3ia_cvc4_pred_all_consts_ ? 2 : 0,
                        options_.ic3ia_cvc4_pred_all_sorts_);
  cvc4a::Grammar g_with_values =
      cvc4_make_grammar(cvc4_solver,
                        cvc4_boundvars,
                        &gs,
                        cvc4_max_terms,
                        options_.ic3ia_cvc4_pred_all_consts_ ? 2 : 1,
                        options_.ic3ia_cvc4_pred_all_sorts_);

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

  logger.log(1, "Looking for {} predicate(s) with CVC4 SyGuS", num_preds);
  try {
    res = cvc4_solver.checkSynth().isUnsat();
    logger.log(1, "Successfully found {} predicate(s)", num_preds);
  }
  catch (cvc4a::CVC4ApiException & e) {
    logger.log(1, "Caught exception from CVC4: {}", e.what());
    return false;
  }

  // for debugging:
  if (options_.verbosity_ > 0) {
    cvc4_solver.printSynthSolution(std::cout);
  }

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
