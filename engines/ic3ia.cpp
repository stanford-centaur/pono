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

// includes for the CVC5 pred experiment
#include "smt-switch/cvc5_solver.h"
#include "smt-switch/cvc5_sort.h"
#include "smt-switch/cvc5_term.h"
#include "smt-switch/identity_walker.h"
#include "smt-switch/utils.h"

#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/term_analysis.h"

using namespace smt;
using namespace std;

namespace pono {

using CVC5SortSet = std::unordered_set<cvc5::Sort>;
using CVC5TermVec = std::vector<cvc5::Term>;

const unordered_set<cvc5::Kind> bv_ops({ cvc5::Kind::EQUAL,
                                          cvc5::Kind::BITVECTOR_CONCAT,
                                          cvc5::Kind::BITVECTOR_EXTRACT,
                                          cvc5::Kind::BITVECTOR_NOT,
                                          cvc5::Kind::BITVECTOR_NEG,
                                          cvc5::Kind::BITVECTOR_AND,
                                          cvc5::Kind::BITVECTOR_OR,
                                          cvc5::Kind::BITVECTOR_XOR,
                                          cvc5::Kind::BITVECTOR_NAND,
                                          cvc5::Kind::BITVECTOR_NOR,
                                          cvc5::Kind::BITVECTOR_XNOR,
                                          cvc5::Kind::BITVECTOR_COMP,
                                          cvc5::Kind::BITVECTOR_ADD,
                                          cvc5::Kind::BITVECTOR_SUB,
                                          cvc5::Kind::BITVECTOR_MULT,
                                          cvc5::Kind::BITVECTOR_UDIV,
                                          cvc5::Kind::BITVECTOR_SDIV,
                                          cvc5::Kind::BITVECTOR_UREM,
                                          cvc5::Kind::BITVECTOR_SREM,
                                          cvc5::Kind::BITVECTOR_SMOD,
                                          cvc5::Kind::BITVECTOR_SHL,
                                          cvc5::Kind::BITVECTOR_ASHR,
                                          cvc5::Kind::BITVECTOR_LSHR,
                                          cvc5::Kind::BITVECTOR_ULT,
                                          cvc5::Kind::BITVECTOR_ULE,
                                          cvc5::Kind::BITVECTOR_UGT,
                                          cvc5::Kind::BITVECTOR_UGE,
                                          cvc5::Kind::BITVECTOR_SLT,
                                          cvc5::Kind::BITVECTOR_SLE,
                                          cvc5::Kind::BITVECTOR_SGT,
                                          cvc5::Kind::BITVECTOR_SGE,
                                          cvc5::Kind::BITVECTOR_ZERO_EXTEND,
                                          cvc5::Kind::BITVECTOR_SIGN_EXTEND,
                                          cvc5::Kind::BITVECTOR_REPEAT,
                                          cvc5::Kind::BITVECTOR_ROTATE_LEFT,
                                          cvc5::Kind::BITVECTOR_ROTATE_RIGHT });

const unordered_set<cvc5::Kind> relational_ops({
    cvc5::Kind::EQUAL,
    cvc5::Kind::DISTINCT,
    cvc5::Kind::LT,
    cvc5::Kind::LEQ,
    cvc5::Kind::GT,
    cvc5::Kind::GEQ,
    cvc5::Kind::BITVECTOR_ULT,
    cvc5::Kind::BITVECTOR_ULE,
    cvc5::Kind::BITVECTOR_UGT,
    cvc5::Kind::BITVECTOR_UGE,
    cvc5::Kind::BITVECTOR_SLT,
    cvc5::Kind::BITVECTOR_SLE,
    cvc5::Kind::BITVECTOR_SGT,
    cvc5::Kind::BITVECTOR_SGE,
});

const unordered_set<cvc5::Kind> multisort_ops({ cvc5::Kind::BITVECTOR_EXTRACT,
                                                 cvc5::Kind::BITVECTOR_CONCAT,
                                                 cvc5::Kind::BITVECTOR_ZERO_EXTEND,
                                                 cvc5::Kind::BITVECTOR_SIGN_EXTEND,
                                                 cvc5::Kind::BITVECTOR_COMP });

const unordered_set<cvc5::Kind> unary_ops({ cvc5::Kind::BITVECTOR_NEG,
                                             cvc5::Kind::BITVECTOR_NOT,
                                             cvc5::Kind::BITVECTOR_EXTRACT,
                                             cvc5::Kind::BITVECTOR_ZERO_EXTEND,
                                             cvc5::Kind::BITVECTOR_SIGN_EXTEND,
                                             cvc5::Kind::NEG });

const unordered_set<cvc5::Kind> bool_ops(
    { cvc5::Kind::AND, cvc5::Kind::OR, cvc5::Kind::XOR, cvc5::Kind::NOT, cvc5::Kind::IMPLIES, cvc5::Kind::ITE });

// Helpers for CVC5 SyGuS Predicate Search
// should eventually be moved elsewhere

bool cvc5_term_is_value(const cvc5::Term & term)
{
  cvc5::Kind k = term.getKind();
  return ((k == cvc5::Kind::CONST_BOOLEAN) || (k == cvc5::Kind::CONST_BITVECTOR)
          || (k == cvc5::Kind::CONST_RATIONAL) || (k == cvc5::Kind::CONST_FLOATINGPOINT)
          || (k == cvc5::Kind::CONST_ROUNDINGMODE) || (k == cvc5::Kind::CONST_STRING)
          || (k == cvc5::Kind::CONST_ARRAY));
}

using CVC5OpSignatures =
    unordered_map<cvc5::Op,
                  unordered_set<cvc5::Sort>>;
using CVC5ValueMap =
    unordered_map<cvc5::Sort,
                  unordered_set<cvc5::Term>>;

/** \class CVC5GrammarSeed
 *  \brief A class for seeding a SyGuS grammar with terms
 *
 *  Used to store ops and values used in an abstract trace
 *  More specifically, it keeps track not only of the operators,
 *  but also which sorts they're applied to
 */
class CVC5GrammarSeed
{
 public:
  CVC5GrammarSeed(cvc5::Solver & solver) : solver_(solver), num_values_(0) {}

  void scan(cvc5::Term term)
  {
    CVC5TermVec to_visit({ term });
    unordered_set<cvc5::Term> visited;
    unordered_set<cvc5::Term> contains_free_vars;
    cvc5::Term t;
    while (!to_visit.empty()) {
      t = to_visit.back();
      cvc5::Sort sort = t.getSort();
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
        // NOTE: can't just use the CVC5 equivalent of is_value
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

        if (cvc5_term_is_value(t)) {
          value_map_[sort].insert(t);
          num_values_++;
        } else if (t.hasOp()) {
          cvc5::Op op = t.getOp();
          // easiest way to store signature is as a function sort
          vector<cvc5::Sort> sort_vec;
          for (const auto & tt : t) {
            sort_vec.push_back(tt.getSort());
          }

          cvc5::Sort ufsort = solver_.mkFunctionSort(sort_vec, sort);
          assert(!op.isNull());
          op_map_[op].insert(ufsort);
        } else if (t.getKind() == cvc5::Kind::CONSTANT) {
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
  const CVC5OpSignatures & get_op_map() const { return op_map_; }

  /** Getter for value map
   *  @return map from sorts to a set of values used of that sort
   */
  const CVC5ValueMap & get_value_map() const { return value_map_; }

  /** Setter for value map
   *  Allows manually setting values in the GrammarSeed
   *  @param value_map the new value map
   */
  void set_value_map(const CVC5ValueMap & value_map) { value_map_ = value_map; }

  const unordered_set<cvc5::Sort> & get_all_sorts()
      const
  {
    return all_sorts_;
  }

  size_t num_values() const { return num_values_; }

 protected:
  cvc5::Solver & solver_;
  CVC5OpSignatures op_map_;
  CVC5ValueMap value_map_;
  unordered_set<cvc5::Sort> all_sorts_;
  size_t num_values_;
};

// helper class for generating grammar for CVC5 SyGuS
cvc5::Grammar cvc5_make_grammar(
    cvc5::Solver & cvc5_solver,
    const CVC5TermVec & cvc5_boundvars,
    const CVC5GrammarSeed * gs,
    unordered_map<cvc5::Sort,
                  unordered_set<cvc5::Term>> & cvc5_max_terms,
    // 0 - no values
    // 1 - values from CVC5GrammarSeed
    // 2 - all values
    size_t values,
    bool all_sorts)
{
  // sorts and their terminal constructors (start constructors)
  cvc5::Sort boolean = cvc5_solver.getBooleanSort();
  cvc5::Term start_bool = cvc5_solver.mkVar(boolean, "Start");
  vector<cvc5::Term> start_terms;

  unordered_map<cvc5::Sort, cvc5::Term> sort2start(
      { { boolean, start_bool } });
  for (const auto & cvc5_bv : cvc5_boundvars)
  {
    cvc5::Sort sort = cvc5_bv.getSort();
    // TODO: remove this limitation
    // TODO-priority: remove bit-vector
    if (!sort.isBitVector() && !sort.isBoolean() && !sort.isReal()
        && !sort.isInteger()) {
      cout << "Skipping unsupported sort " << sort << endl;
      continue;
    }

    // AI: keep this...
    if (sort2start.find(sort) == sort2start.end())
    {
      cvc5::Term start_term = cvc5_solver.mkVar(sort, sort.toString() + "_start");
      sort2start[sort] = start_term;
      start_terms.push_back(start_term);
    }
  }

  unordered_map<cvc5::Term, vector<cvc5::Term>>
      constructs;

  // collect variables by sort and add to constructs
  unordered_map<cvc5::Sort, vector<cvc5::Term>>
      sort2boundvars;
  for (const auto & bvar : cvc5_boundvars) {
    cvc5::Sort sort = bvar.getSort();
    sort2boundvars[sort].push_back(bvar);
    constructs[sort2start.at(sort)].push_back(bvar);
  }

  // TODO: remove this grammar seed
  if (false && gs) {
    const CVC5OpSignatures & ops_map = gs->get_op_map();
    const CVC5ValueMap & values_map = gs->get_value_map();

    if (all_sorts) {
      // add all sorts to start terms
      for (const auto & sort : gs->get_all_sorts()) {
        // treat booleans specially, only included as result of relational
        // operations
        if (sort2start.find(sort) == sort2start.end() && !sort.isBoolean()) {
          cvc5::Term start_term =
              cvc5_solver.mkVar(sort, sort.toString() + "_start");
          sort2start[sort] = start_term;
          start_terms.push_back(start_term);
        }
      }
    }

    // separate operators
    unordered_map<SortKind, unordered_set<cvc5::Kind>> reg_ops;
    unordered_map<SortKind, unordered_set<cvc5::Kind>> rel_ops;
    unordered_map<SortKind, unordered_set<cvc5::Op>>
        ms_ops;

    assert(ops_map.size());
    for (const auto & opelem : ops_map) {
      cvc5::Op op = opelem.first;
      cvc5::Kind po = op.getKind();
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
      cvc5::Sort sort = s.getSort();
      SortKind sk;
      if (sort.isBitVector()) {
        sk = BV;
      } else if (sort.isReal() || sort.isInteger()) {
        // using real for both integers and reals
        sk = REAL;
        logger.log(0, "Adding arithmetic ADD\n");
        constructs[s].push_back(cvc5_solver.mkTerm(cvc5::Kind::ADD, {s, s}));
        //        constructs[s].push_back(cvc5_solver.mkTerm(cvc5::MULT, {s, s}));
      } else if (sort.isBoolean()) {
        // only looking for predicates
        // nothing to do with a boolean
        // shouldn't be included in start terms
        assert(false);
      } else {
        // TODO handle INT/REAL (note using REAL for both)
        logger.log(0,
                   "WARNING IC3IA CVC5 predicates unhandled sort: {}",
                   s.getSort().toString());
      }

      // regular
      for (const auto & po : reg_ops[sk]) {
        if (unary_ops.find(po) == unary_ops.end()) {
          constructs[s].push_back(cvc5_solver.mkTerm(po, {s, s}));
        } else {
          constructs[s].push_back(cvc5_solver.mkTerm(po, {s}));
        }
      }

      // relational
      for (const auto & po : rel_ops[sk]) {
        constructs[start_bool].push_back(cvc5_solver.mkTerm(po, {s, s}));
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
          vector<cvc5::Term> args;
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

          cvc5::Sort codomain_sort = signature.getFunctionCodomainSort();
          if (sort2start.find(codomain_sort) != sort2start.end()) {
            cvc5::Term return_start_term = sort2start.at(codomain_sort);
            constructs[return_start_term].push_back(
                cvc5_solver.mkTerm(op, args));
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
      //if ( s is of sort BV) {
      cvc5::Term equals = cvc5_solver.mkTerm(cvc5::Kind::EQUAL, {s, s});
      cvc5::Term bvugt = cvc5_solver.mkTerm(cvc5::Kind::BITVECTOR_UGT, {s, s});
      constructs[start_bool] = { bvugt, equals };
      //} else if (s is of sort int/real) {
      //cvc5::Term equals = cvc5_solver.mkTerm(cvc5::EQUAL, {s, s});
      //cvc5::Term bvugt = cvc5_solver.mkTerm(cvc5::BITVECTOR_UGT, {s, s});
      //}
    }

    // include bv operations in the grammar
    for (auto s : start_terms) {
      cvc5::Term zero = cvc5_solver.mkBitVector(s.getSort().getBitVectorSize(), 0);
      cvc5::Term one = cvc5_solver.mkBitVector(s.getSort().getBitVectorSize(), 1);
      cvc5::Term min_signed = cvc5_solver.mkBitVector(s.getSort().getBitVectorSize(), pow(2,s.getSort().getBitVectorSize() - 1));
      cvc5::Term bvadd = cvc5_solver.mkTerm(cvc5::Kind::BITVECTOR_ADD, {s, s});
      cvc5::Term bvmul = cvc5_solver.mkTerm(cvc5::Kind::BITVECTOR_MULT, {s, s});
      cvc5::Term bvand = cvc5_solver.mkTerm(cvc5::Kind::BITVECTOR_AND, {s, s});
      cvc5::Term bvcomp = cvc5_solver.mkTerm(cvc5::Kind::BITVECTOR_COMP, {s, s});
      cvc5::Term bvor = cvc5_solver.mkTerm(cvc5::Kind::BITVECTOR_OR, {s, s});
      cvc5::Term bvnot = cvc5_solver.mkTerm(cvc5::Kind::BITVECTOR_NOT, {s});
      cvc5::Term bvneg = cvc5_solver.mkTerm(cvc5::Kind::BITVECTOR_NEG, {s});
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

  // TODO: remove this
  // add max terms
  for (const auto & start_term : start_terms) {
    cvc5::Sort sort = start_term.getSort();
    for (const auto & t : cvc5_max_terms[sort]) {
      constructs[start_term].push_back(t);
    }
  }

  for (const auto & start_term : start_terms) {
    cvc5::Sort sort = start_term.getSort();
    if (sort2boundvars.find(sort) == sort2boundvars.end()) {
      // there should be some kind of start (e.g. value or constant)
      // for every start term
      // it's possible there will be no constants if the only terms
      // of a given sort contain input variables
      if (!sort.isBitVector()) {
        throw PonoException("Unhandled empty start term of sort: "
                            + sort.toString());
      }

      cvc5::Term zero = cvc5_solver.mkBitVector(sort.getBitVectorSize(), 0);
      cvc5::Term one = cvc5_solver.mkBitVector(sort.getBitVectorSize(), 1);
      constructs[start_term].push_back(zero);
      constructs[start_term].push_back(one);
    }
    // add it for int/real
  }

  // construct the grammar
  vector<cvc5::Term> starts({ start_bool });
  starts.reserve(start_terms.size() + 1);
  starts.insert(starts.end(), start_terms.begin(), start_terms.end());
  cvc5::Grammar g = cvc5_solver.mkGrammar(cvc5_boundvars, starts);

  for (const auto & elem : constructs) {
    const cvc5::Term & start_term = elem.first;
    const vector<cvc5::Term> & rules = elem.second;
    assert(rules.size());

    // add rules for this nonterminal start term
    g.addRules(start_term, rules);

    if (values == 2) {
      g.addAnyConstant(start_term);
    }
  }

  return g;
}

// helper class for translating back from CVC5
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
      ia_(conc_ts_, ts_, opt.ic3ia_cvc5_pred_ ? abs_unroller_ : unroller_),
      // only mathsat interpolator supported
      interpolator_(create_interpolating_solver_for(
          SolverEnum::MSAT_INTERPOLATOR, Engine::IC3IA_ENGINE)),
      to_interpolator_(interpolator_),
      to_solver_(solver_),
      longest_cex_length_(0),
      abs_unroller_(ts_)
{
  // since we passed a fresh RelationalTransitionSystem as the main TS
  // need to point orig_ts_ to the right place
  orig_ts_ = ts;
  engine_ = Engine::IC3IA_ENGINE;
  approx_pregen_ = true;
}

void IC3IA::add_important_var(Term v)
{
  if (!options_.ic3ia_track_important_vars_) {
    return;
  }

  // have to consider that original solver
  // might not be the same as the prover solver
  if (solver_ != orig_ts_.solver()) {
    v = to_prover_solver_.transfer_term(v);
  }
  logger.log(1, "Adding important variable: {}", v);
  ia_.add_important_var(v);
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

  // gather values
  get_leaves(conc_ts_.init(), ts_values_);
  get_leaves(conc_ts_.trans(), ts_values_);
  get_leaves(bad_, ts_values_);
  // remove all non-ts_values_ (leaves include symbols)
  // TODO create helper function for this
  for (auto it = ts_values_.begin(); it != ts_values_.end();) {
    if (!(*it)->is_value()) {
      it = ts_values_.erase(it);
    } else {
      it++;
    }
  }

  // gather max-terms (for ic3ia-cvc5-pred)
  if (options_.ic3ia_cvc5_pred_maxterms_) {
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

  // HACK added to experiment with CVC5 SyGuS for finding predicates
  if (options_.ic3ia_cvc5_pred_) {
    // first check if the cex trace is spurious
    // this is a bit hacky for now -- should refactor later so that
    // the unrolling is shared with the interpolator approach in the else branch
    // but don't want to disentagle it from interpolation right now
    assert(solver_context_ == 0);
    push_solver_context();

    // NOTE only use abs_unroller in ic3ia-cvc5-pred mode
    // need it to unroll the ts_ in cvc5_find_preds
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
      cvc5_find_preds(cex_, preds);
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

    // new predicates
    for (auto const&p : preds) {
      if (predset_.find(p) == predset_.end()) {
        // unseen predicate
        fresh_preds.push_back(p);
      }
    }
  }

  if (!fresh_preds.size()) {
    logger.log(1, "IC3IA: refinement failed couldn't find any new predicates");
    return RefineResult::REFINE_FAIL;
  }

  if (options_.random_seed_ > 0) {
    shuffle(fresh_preds.begin(),
            fresh_preds.end(),
            default_random_engine(options_.random_seed_));
  }

  // reduce new predicates
  TermVec red_preds;
  if (options_.ic3ia_reduce_preds_
      && ia_.reduce_predicates(cex_, fresh_preds, red_preds)) {
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
    // if enabled should only fail if removed all predicates
    // this can happen when there are uninterpreted functions
    // the unrolling can force incompatible UF interpretations
    // but IC3 (which doesn't unroll) still needs the predicates
    // in this case, just use all the fresh predicates
    assert(!options_.ic3ia_reduce_preds_ || red_preds.size() == 0);
    logger.log(2, "reduce predicates FAILED");
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

bool IC3IA::cvc5_find_preds(const TermVec & cex, UnorderedTermSet & out_preds)
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
      // by traversing cvc5_formula (at smt-switch level)
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
  size_t num_preds = 3;
  while (!found_preds) {
    found_preds = cvc5_synthesize_preds(
        abs_trace, statevars, var_args, free_vars, num_preds, out_preds);
    num_preds++;
  }

  return found_preds;
}

bool IC3IA::cvc5_synthesize_preds(
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
  // hacked in to evaluate CVC5
  // if done for real, should be sure to do this OR the interpolator, not both
  smt::SmtSolver cvc5_ = create_solver(smt::SolverEnum::CVC5);
  smt::TermTranslator to_cvc5_(cvc5_);
  smt::TermTranslator from_cvc5_(solver_);

  // populate from_cvc5_ cache
  UnorderedTermMap & from_cvc5_cache = from_cvc5_.get_cache();
  for (auto sv : conc_ts_.statevars()) {
    from_cvc5_cache[to_cvc5_.transfer_term(sv)] = sv;
  }

  // get the underlying CVC5 objects
  cvc5::Term cvc5_formula =
      static_pointer_cast<Cvc5Term>(to_cvc5_.transfer_term(abs_trace, BOOL))
          ->get_cvc5_term();

  unordered_set<cvc5::Term> cvc5_free_vars;
  for (auto fv : free_vars) {
    cvc5_free_vars.insert(
        static_pointer_cast<Cvc5Term>(to_cvc5_.transfer_term(fv))
            ->get_cvc5_term());
  }

  // Retrieve the underlying fresh cvc5 solver
  cvc5::Solver & cvc5_solver =
      static_pointer_cast<Cvc5Solver>(cvc5_)->get_cvc5_solver();

  // set necessary options for sygus
  cvc5_solver.setOption("sygus", "true");
  cvc5_solver.setOption("lang", "sygus2");
  // we can set the incremental mode to true
  cvc5_solver.setOption("incremental", "false");

  if (options_.ic3ia_cvc5_pred_size_)
  {
    int pred_size = options_.ic3ia_cvc5_pred_size_;
    pred_size += (num_preds-1)/3; // increase the size periodically
    // sygus-abort-size: number of operators (+-1)
    cvc5_solver.setOption("sygus-abort-size", std::to_string(pred_size));
  }
  // TODO: what is the new option for sygus enum mode?
  cvc5_solver.setOption("sygus-enum", "fast");

  // create bound variables to use in the synthesized function
  vector<cvc5::Term> cvc5_statevars;
  cvc5_statevars.reserve(statevars.size());
  vector<cvc5::Term> cvc5_boundvars;
  cvc5_boundvars.reserve(statevars.size());
  Term transferred_sv;
  for (auto sv : statevars) {
    transferred_sv = to_cvc5_.transfer_term(sv);
    cvc5::Term cvc5_sv =
        static_pointer_cast<Cvc5Term>(transferred_sv)->get_cvc5_term();
    cvc5_statevars.push_back(cvc5_sv);
    // mkVar is for Sygus
    cvc5::Term cvc5_bv =
        cvc5_solver.mkVar(cvc5_sv.getSort(), cvc5_sv.toString() + "_var");
    cvc5_boundvars.push_back(cvc5_bv);
  }

  // TODO: remove this max-term thing
  unordered_map<cvc5::Sort,
                unordered_set<cvc5::Term>>
      cvc5_max_terms;
  for (const auto & t : max_terms_) {
    cvc5::Term cvc5_t =
        static_pointer_cast<Cvc5Term>(to_cvc5_.transfer_term(t))
            ->get_cvc5_term();
    // need to substitute with bound vars
    cvc5_max_terms[cvc5_t.getSort()].insert(
        cvc5_t.substitute(cvc5_statevars, cvc5_boundvars));
  }

  CVC5GrammarSeed gs(cvc5_solver);
  Term transferred_trace = to_cvc5_.transfer_term(abs_trace, BOOL);
  gs.scan(static_pointer_cast<Cvc5Term>(transferred_trace)->get_cvc5_term());

  logger.log(1, "Number of Values in Trace: {}", gs.num_values());

  // Grammar construction
  cvc5::Grammar g_with_trace_values =
      cvc5_make_grammar(cvc5_solver,
                        cvc5_boundvars,
                        &gs,
                        cvc5_max_terms,
                        options_.ic3ia_cvc5_pred_all_consts_ ? 2 : 1,
                        options_.ic3ia_cvc5_pred_all_sorts_);

  // prioritize using only values in ts first (not the trace
  // which could have new values due to top-level substitution)
  CVC5ValueMap cvc5_value_map;
  for (const auto & val : ts_values_) {
    cvc5::Term cvc5_val =
        static_pointer_cast<Cvc5Term>(to_cvc5_.transfer_term(val))
            ->get_cvc5_term();
    cvc5_value_map[cvc5_val.getSort()].insert(cvc5_val);
  }
  gs.set_value_map(cvc5_value_map);
  cvc5::Grammar g =
      cvc5_make_grammar(cvc5_solver,
                        cvc5_boundvars,
                        &gs,
                        cvc5_max_terms,
                        options_.ic3ia_cvc5_pred_all_consts_ ? 2 : 1,
                        options_.ic3ia_cvc5_pred_all_sorts_);

  vector<cvc5::Term> pred_vec;
  for (size_t n = 0; n < num_preds; ++n) {
    // Create the predicate to search for. Use the grammar
    string pred_name = "P_" + std::to_string(n);
    // cvc5::Term pred = 
    //   cvc5_solver.synthFun(pred_name, cvc5_boundvars, cvc5_solver.getBooleanSort(), g);
    cvc5::Term pred;
    switch (n % 3) {
    case 0:
      pred = cvc5_solver.synthFun(pred_name, cvc5_boundvars,
                                  cvc5_solver.getBooleanSort(), g);
      break;
    case 1:
      pred = cvc5_solver.synthFun(pred_name,
                                  cvc5_boundvars,
                                  cvc5_solver.getBooleanSort(),
                                  g_with_trace_values);
      break;
    default:
      pred = cvc5_solver.synthFun(pred_name, cvc5_boundvars,
                                  cvc5_solver.getBooleanSort());
      break;
    };
    pred_vec.push_back(pred);

    // add the implicit predicate abstraction constraints
    // e.g. P(x^@0) <-> P(x@1) /\ P(x^@1) <-> P(x@2) /\ ...

    vector<cvc5::Term> cvc5_next_var_args;
    vector<cvc5::Term> cvc5_abs_var_args;

    for (const auto & var_pair : unrolled_var_args) {
      const TermVec & unrolled_next_vars = var_pair.first;
      const TermVec & unrolled_abs_vars = var_pair.second;

      // CVC5 takes the function as the first child, so insert the
      // pred function first
      cvc5_next_var_args.clear();
      cvc5_next_var_args.push_back(pred);

      cvc5_abs_var_args.clear();
      cvc5_abs_var_args.push_back(pred);

      assert(unrolled_next_vars.size() == unrolled_abs_vars.size());
      Term nv, abs_nv;
      for (size_t i = 0; i < unrolled_next_vars.size(); ++i) {
        nv = to_cvc5_.transfer_term(unrolled_next_vars[i]);
        abs_nv = to_cvc5_.transfer_term(unrolled_abs_vars[i]);

        cvc5_next_var_args.push_back(
            static_pointer_cast<Cvc5Term>(nv)->get_cvc5_term());
        cvc5_abs_var_args.push_back(
            static_pointer_cast<Cvc5Term>(abs_nv)->get_cvc5_term());
      }

      cvc5::Term pred_app_vars =
          cvc5_solver.mkTerm(cvc5::Kind::APPLY_UF, cvc5_next_var_args);
      cvc5::Term pred_app_abs_vars =
          cvc5_solver.mkTerm(cvc5::Kind::APPLY_UF, cvc5_abs_var_args);

      cvc5_formula = cvc5_solver.mkTerm(cvc5::Kind::AND,
                                        {cvc5_formula,
                                         cvc5_solver.mkTerm(cvc5::Kind::EQUAL, {pred_app_vars, pred_app_abs_vars})});
    }
  }

  cvc5::Term constraint = cvc5_solver.mkTerm(cvc5::Kind::NOT, {cvc5_formula});

  // use sygus variables rather than ordinary variables.
  std::map<cvc5::Term, cvc5::Term> old_to_new;
  std::vector<cvc5::Term> originals(cvc5_free_vars.begin(),
                                     cvc5_free_vars.end());
  for (cvc5::Term old_var : originals) {
    cvc5::Term new_var =
      cvc5_solver.declareSygusVar(old_var.toString() + "_sy", old_var.getSort());
    old_to_new[old_var] = new_var;
  }
  std::vector<cvc5::Term> news;
  for (cvc5::Term old_var : originals) {
    assert(old_to_new.find(old_var) != old_to_new.end());
    news.push_back(old_to_new[old_var]);
  }

  cvc5::Term sygus_constraint = constraint.substitute(originals, news);
  cvc5_solver.addSygusConstraint(sygus_constraint);

  logger.log(1, "Looking for {} predicate(s) with CVC5 SyGuS", num_preds);
  try {
    res = cvc5_solver.checkSynth().hasSolution();
    logger.log(1, "Successfully found {} predicate(s)", num_preds);
  }
  catch (cvc5::CVC5ApiException & e) {
    logger.log(1, "Caught exception from CVC5: {}", e.what());
    return false;
  }

  // for debugging:
  if (options_.verbosity_ > 0) {
    // TODO: fix printout on std::cout
    //cvc5_solver.getSynthSolution();
  }

  for (auto pred : pred_vec) {
    cvc5::Term pred_solution = cvc5_solver.getSynthSolution(pred);

    // instead of applying predicate and simplifying to get rid of the lambda
    // decided to just do substitution
    // this keeps the structure that sygus came up with instead of rewriting it
    // noticed on one example that arithmetic bv operators were replaced with
    // bitwise operators

    vector<cvc5::Term> pred_solution_children(pred_solution.begin(),
                                               pred_solution.end());
    assert(pred_solution_children.size()
           == 2);  // expecting a bound_var_list and a term
    cvc5::Term pred_solution_statevars =
      pred_solution_children[1].substitute(cvc5_boundvars, cvc5_statevars);

    Term cvc5_learned_pred = make_shared<Cvc5Term>(pred_solution_statevars);
    // Makai: put this back in if it starts failing on term translation
    //        not all backend solvers support n-ary arguments
    // Binarizer binarizer(cvc5_);
    // cvc5_learned_pred = binarizer.process(cvc5_learned_pred);

    Term learned_pred = from_cvc5_.transfer_term(cvc5_learned_pred, BOOL);
    assert(learned_pred);

    // NOTE in future might look for more than one predicate at a time
    //out_preds.insert(learned_pred);
    get_predicates(solver_, learned_pred, out_preds, false, false, true);
  }

  return res;
}

}  // namespace pono
