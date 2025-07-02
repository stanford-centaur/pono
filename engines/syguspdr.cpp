/*********************                                                  */
/*! \file syguspdr.cpp
** \verbatim
** Top contributors (to current version):
**   Hongce Zhang
** This file is part of the pono project.
** Copyright (c) 2020 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Syntax-Guided Synthesis in IC3/PDR implementation
**
**
**/

#include "engines/syguspdr.h"

#include "modifiers/mod_ts_prop.h"
#include "utils/container_shortcut.h"
#include "utils/logger.h"
#include "utils/sygus_predicate_constructor.h"
#include "utils/term_walkers.h"

using namespace smt;

namespace pono {

// #define DEBUG
#ifdef DEBUG
#define D(...) logger.log(__VA_ARGS__)
#define INFO(...) D(0, __VA_ARGS__)
#else
#define D(...) \
  do {         \
  } while (0)
#define INFO(...) logger.log(3, __VA_ARGS__)
#endif

// ----------------------------------------------------------------
// Helper functions
// ----------------------------------------------------------------

unsigned SygusPdr::GetScore(const Term & t)
{
  const auto & score_map = term_score_walker_.GetScoreMap();
  auto pos = score_map.find(t);
  if (pos != score_map.end()) return pos->second.score;
  term_score_walker_.WalkBFS(t);
  pos = score_map.find(t);
  assert(pos != score_map.end());
  return pos->second.score;
}

// call this before op abstraction
// op abstraction will change ts_.inputvar?
// and promote again is needed, but we want to
// keep these states in the model
void SygusPdr::build_ts_related_info()
{
  // the input vars and the prime to next function
  const auto & all_state_vars = ts_.statevars();
  for (const auto & sv : all_state_vars) {
    const auto & s_updates = ts_.state_updates();
    if (!IN(sv, s_updates))
      no_next_vars_.insert(sv);
    else
      nxt_state_updates_.emplace(ts_.next(sv), s_updates.at(sv));
  }
}

smt::Term SygusPdr::next_curr_replace(const smt::Term & in) const
{
  return ts_.solver()->substitute(in, nxt_state_updates_);
}

bool SygusPdr::test_ts_has_op(const std::unordered_set<PrimOp> & prim_ops) const
{
  UnorderedTermSet term_op_out;
  TermOpCollector op_collector(ts_.solver());
  for (const auto & s_update : ts_.state_updates()) {
    op_collector.find_matching_terms(s_update.second, prim_ops, term_op_out);
  }
  return !term_op_out.empty();
}

// ----------------------------------------------------------------
// Constructor & Destructor
// ----------------------------------------------------------------
SygusPdr::SygusPdr(const SafetyProperty & p,
                   const TransitionSystem & ts,
                   const SmtSolver & s,
                   PonoOptions opt)
    : super(p, ts, s, opt),
      partial_model_getter(solver_),
      has_assumptions(true)  // most conservative way
{
  solver_->set_opt("produce-unsat-assumptions", "true");

  // we need to have the reset-assertion capability
  if (solver_->get_solver_enum() == SolverEnum::BTOR)
    solver_->set_opt("base-context-1", "true");
  // actually the above is already done in
}

// destruct -- free the memory
SygusPdr::~SygusPdr()
{
  for (const auto & ic3f_mptr : model2cube_) {
    if (ic3f_mptr.second) delete ic3f_mptr.second;
  }

  //  for (const auto & partial2full : to_full_model_map_) {
  //    if (partial2full.second)
  //      delete partial2full.second; // delete the full models
  //  }
}  // SygusPdr::~SygusPdr

// ----------------------------------------------------------------
// IC3 override
// ----------------------------------------------------------------

// init
void SygusPdr::initialize()
{
  if (!options_.ic3_pregen_)
    throw PonoException("SyGuS PDR requires IC3 predecessor generalization");
  if (!options_.ic3_indgen_)
    throw PonoException("SyGuS PDR requires IC3 inductive generalization");
  if (options_.assume_prop_)
    throw PonoException(
        "SyGuS PDR requires not to assume property in the transition "
        "relations");
  options_.ic3_unsatcore_gen_ = false;

  if (initialized_) return;

  ts_ = orig_ts_;  // maybe call promote inputvars implicitly here
  // I really need the prime variable for inputs
  // otherwise the corner cases are hard to handle...

  super::initialize();  // I don't need the trans->prop thing

  bad_next_ = ts_.next(bad_);  // bad is only available after parent's init

  if (options_.sygus_use_operator_abstraction_) {
    if (options_.sygus_use_operator_abstraction_ == 2)
      op_abstractor_ = std::make_unique<OpUfAbstractor>(
          orig_ts_,
          ts_,
          std::unordered_set<PrimOp>(
              { BVMul, BVUdiv, BVSdiv, BVSmod, BVSrem, BVUrem }));
    else {
      op_abstractor_ = std::make_unique<OpInpAbstractor>(
          orig_ts_,
          ts_,
          std::unordered_set<PrimOp>(
              { BVMul, BVUdiv, BVSdiv, BVSmod, BVSrem, BVUrem }),
          bad_,
          0);

      ts_ = promote_inputvars(ts_);
    }
    // op_abstractor_ = std::make_unique<OpUfAbstractor>(
    //     orig_ts_,
    //     ts_,
    //     std::unordered_set<PrimOp>(
    //         { BVMul, BVUdiv, BVSdiv, BVSmod, BVSrem, BVUrem }),
    //     bad_,
    //     4);
    if (options_.sygus_term_mode_ == SyGuSTermMode::TERM_MODE_AUTO)
      options_.sygus_term_mode_ = op_abstractor_->has_abstracted()
                                      ? (SyGuSTermMode::SPLIT_FROM_DESIGN)
                                      : (SyGuSTermMode::FROM_DESIGN_LEARN_EXT);
    // we need to reset trans function in the base class
    reset_solver();
  } else {
    if (options_.sygus_term_mode_ == SyGuSTermMode::TERM_MODE_AUTO) {
      options_.sygus_term_mode_ =
          (test_ts_has_op({ BVMul, BVUdiv, BVSdiv, BVSmod, BVSrem, BVUrem })
               ? SyGuSTermMode::SPLIT_FROM_DESIGN
               : (test_ts_has_op({ BVAdd, BVSub })
                      ? SyGuSTermMode::FROM_DESIGN_LEARN_EXT
                      : (test_ts_has_op({ BVUle, BVUlt })
                             ? SyGuSTermMode::VAR_C_EXT
                             : SyGuSTermMode::VAR_C_EXT)));
    }
  }

  if (options_.sygus_term_mode_ == SyGuSTermMode::TERM_MODE_AUTO)
    options_.sygus_term_mode_ = SyGuSTermMode::FROM_DESIGN_LEARN_EXT;

  build_ts_related_info();

  // has_assumption -- on the original one
  has_assumptions = false;
  assert(!nxt_state_updates_.empty());
  for (const auto & c_initnext : ts_.constraints()) {
    // if (!c_initnext.second)
    //  continue; // should not matter
    has_assumptions = true;
    assert(ts_.no_next(c_initnext.first));
    // if (no_next) {
    constraints_curr_var_.push_back(c_initnext.first);
    // translate input_var to next input_var
    // but the state var ...
    // we will get to next anyway
    constraints_curr_var_.push_back(
        next_curr_replace(ts_.next(c_initnext.first)));
    // } // else skip
  }
  // initialize the caches
  // extract the operators
  op_extract_ = std::make_unique<syntax_analysis::OpExtractor>();
  op_extract_->WalkBFS(ts_.init());
  op_extract_->WalkBFS(ts_.trans());
  op_extract_->GetSyntaxConstruct().RemoveConcat();
  op_extract_->GetSyntaxConstruct().RemoveExtract();
  op_extract_->GetSyntaxConstruct().AndOrConvert();
  op_extract_->GetSyntaxConstruct().RemoveUnusedStructure();

  {  // 1. register terms to find exprs
    // 2. extract parent from the same terms
    for (auto && v_nxtexpr_pair : ts_.state_updates()) {
      assert(ts_.no_next(v_nxtexpr_pair.second));
      sygus_term_manager_.RegisterTermsToWalk(v_nxtexpr_pair.second);
      parent_of_terms_.WalkBFS(v_nxtexpr_pair.second);
    }
    sygus_term_manager_.RegisterTermsToWalk(ts_.init());
    parent_of_terms_.WalkBFS(ts_.init());

    for (const auto & c_next_init : ts_.constraints()) {
      if (!c_next_init.second) continue;
      assert(ts_.no_next(c_next_init.first));
      // if (ts_.no_next(c_next_init.first)) {
      sygus_term_manager_.RegisterTermsToWalk(c_next_init.first);
      parent_of_terms_.WalkBFS(c_next_init.first);
      //  }
    }

    // sygus_term_manager_.RegisterTermsToWalk(property_.prop());
    // parent_of_terms_.WalkBFS(property_.prop());
  }

  // cache two lambda functions for sygus enum
  to_next_func_ = [this](const Term & v) -> Term { return this->ts_.next(v); };

  score_func_ = [this](const Term & v) -> unsigned {
    return this->GetScore(v);
  };

  {  // now create TermLearner
    term_learner_.reset(new syntax_analysis::TermLearner(
        to_next_func_,
        score_func_,
        solver_,  // you need to make trans dynamic
        parent_of_terms_));
  }

}  // initialize

void SygusPdr::check_ts() const  // custom_ts_ is not ready at this point
{
  // check if there are arrays or uninterpreted sorts and fail if so
  if (!ts_.is_functional())
    throw PonoException(
        "SyGuS PDR only supports functional transition systems.");
  // check if there are arrays or uninterpreted sorts and fail if so
  for (const auto & vec : { ts_.statevars(), ts_.inputvars() }) {
    for (const auto & st : vec) {
      SortKind sk = st->get_sort()->get_sort_kind();
      if (sk == ARRAY) {
        throw PonoException("SyGuS PDR does not support arrays yet");
      } else if (sk == UNINTERPRETED) {
        throw PonoException(
            "SyGuS PDR does not support uninterpreted sorts yet.");
      }
    }
  }

  if (!ts_.inputvars().empty()) {
    throw PonoException(
        "SyGuS PDR requires promoting input variables to state variables.");
  }
}  // check_ts

bool SygusPdr::ic3formula_check_valid(const IC3Formula & u) const
{
  const Sort & boolsort = solver_->make_sort(BOOL);
  // check that children are literals
  Term pred;
  Op op;
  for (const auto & c : u.children) {
    if (c->get_sort() != boolsort) {
      logger.log(
          3, "ERROR SyGuS PDR IC3Formula contains non-boolean atom: {}", c);
      return false;
    }
  }
  // got through all checks without failing
  return true;
}

IC3Formula SygusPdr::get_model_ic3formula() const
{
  // generalize_predecessor does not use its output at all
  Term conj_full;
  TermVec conjvec;
  for (const auto & v : ts_.statevars()) {
    Term val = solver_->get_value(v);
    auto eq = solver_->make_term(Op(PrimOp::Equal), v, val);
    conjvec.push_back(eq);

    if (conj_full)
      conj_full = solver_->make_term(Op(PrimOp::And), conj_full, eq);
    else
      conj_full = eq;
  }

  for (const auto & v : ts_.inputvars()) {
    Term val = solver_->get_value(v);
    auto eq = solver_->make_term(Op(PrimOp::Equal), v, val);
    conjvec.push_back(eq);

    if (conj_full)
      conj_full = solver_->make_term(Op(PrimOp::And), conj_full, eq);
    else
      conj_full = eq;
  }
  assert(conj_full);
  assert(!conjvec.empty());

  return IC3Formula(conj_full, { conj_full }, false /*not a disjunction*/);
}  // SygusPdr::get_model_ic3formula

IC3Formula SygusPdr::inductive_generalization(size_t i, const IC3Formula & c)
{
  assert(i > 0);
  assert(!c.disjunction);  // we need it to be a model

  // this is to produce the lemma
  // find the model
  auto model_pos = model2cube_.find(c.term);
  assert(model_pos != model2cube_.end());
  syntax_analysis::IC3FormulaModel * post_model = model_pos->second;

  Term Init_prime = ts_.next(ts_.init());
  const Term & cex = c.term;
  // Term not_cex = smart_not(cex);
  Term cex_prime = ts_.next(cex);
  Term Fprev = get_frame_term(i - 1);
  // since we add the property as a lemma, so we should be okay here
  Term T = ts_.trans();
  Term F_T_not_cex = make_and({ Fprev, T });  // , not_cex
  Term base = solver_->make_term(Or, F_T_not_cex, Init_prime);

  IC3Formula pre_formula;
  syntax_analysis::IC3FormulaModel * pre_model = NULL;
  syntax_analysis::IC3FormulaModel * pre_full_model = NULL;
  bool failed_at_init;

  bool insufficient_pred = false;
  do {
    // refresh the frame lemmas
    Fprev = get_frame_term(i - 1);
    T = ts_.trans();
    F_T_not_cex = make_and({ Fprev, T });  // , not_cex
    base = solver_->make_term(Or, F_T_not_cex, Init_prime);

    syntax_analysis::PredConstructor pred_collector_(
        to_next_func_,
        solver_,
        post_model,  // IC3FormulaModel * cex
        setup_cex_info(post_model));

    uint64_t old_term_count = pred_collector_.term_summary();
    {  // step 1 - check if the pred are good if not construct the model
      const auto & pred_nxt_ = pred_collector_.GetAllPredNext();
      push_solver_context();
      disable_all_labels();
      // check ( F[i-1] /\ NOT(c) /\ T ) \/ init' AND pred_next sat

      solver_->assert_formula(base);
      for (const auto & pred_ : pred_nxt_) solver_->assert_formula(pred_);

      Result r = check_sat();
      if (r.is_sat()) {
        // TODO: extract model and do recursive block
        // pre_model will be put into the map and free at class destructor

        insufficient_pred = true;
        failed_at_init =  // if model(init_prime) == 1
            !(syntax_analysis::eval_val(
                  solver_->get_value(Init_prime)->to_string())
              < syntax_analysis::eval_val("#b1"));
        if (pre_full_model) {
          delete pre_full_model;
          pre_full_model = NULL;
        }
        if (!failed_at_init) {
          // the idea of full model is we want to keep the assignment to
          // some inputs that could be thrown away
          std::tie(pre_formula, pre_model, pre_full_model) =
              ExtractPartialAndFullModel(post_model);
          D(3, "Generated MAY-block model {}", pre_model->to_string());
        } else {
          // only pre_model and failed_at_init are used in later
          std::tie(pre_formula, pre_model) = ExtractInitPrimeModel(
              Init_prime);  // this can be only the prime vars
          D(3, "Generated MAY-block model (init) {}", pre_model->to_string());
        }
        // note Here you must give a formula with current variables
      } else
        insufficient_pred = false;
      pop_solver_context();
    }  // end of step 1

    if (insufficient_pred) {
      // extract Proof goal from partial_model
      if (!failed_at_init) {
        D(3, "Try recursive block the above model.");
        if (try_recursive_block_goal(pre_formula, i - 1))
          continue;  // see if we have next that we may need to block
        if (pred_collector_.term_summary() > old_term_count) {
          D(3,
            "though not blocked, we got {} > {} terms",
            pred_collector_.term_summary(),
            old_term_count);
          continue;
        }
      }  // if we failed at init, then we will anyway need to
      // propose new terms
      D(3, "Cannot block MAY-block model {}", pre_model->to_string());
      D(3, "Back to F[{}]->[{}], get new terms.", i - 1, i);

      if (!failed_at_init)
        assert(pre_full_model);
      else
        assert(!pre_full_model);

      bool succ = propose_new_terms(failed_at_init ? pre_model : pre_full_model,
                                    post_model,
                                    F_T_not_cex,
                                    Init_prime,
                                    failed_at_init);
      assert(succ);
      delete pre_full_model;
      pre_full_model = NULL;

      // it has loop inside
      // and it must succeed because eventually we will end with bit-level
      // things
    } else {  // if sufficient
      if (pre_full_model) delete pre_full_model;

      // at this point we have enough preds
      // reduce the preds
      IC3Formula ret =
          options_.smt_solver_ == SolverEnum::BTOR
              ? select_predicates_btor(base, pred_collector_.GetAllPredNext())
              : select_predicates_generic(base,
                                          pred_collector_.GetAllPredNext());

      return ret;
    }
  } while (insufficient_pred);
  assert(false);  // should not be reachable
}  // inductive_generalization

// special optimization for btor
IC3Formula SygusPdr::select_predicates_btor(const Term & base,
                                            const TermVec & preds_nxt)
{
  D(3, "[IterativeReduction] Initial core size {}", preds_nxt.size());
  UnorderedTermSet unsatcore(preds_nxt.begin(), preds_nxt.end());

  push_solver_context();
  disable_all_labels();
  syntax_analysis::reduce_unsat_core_to_fixedpoint(base, unsatcore, solver_);
  pop_solver_context();

  D(3, "[IterativeReduction] End core size {}", unsatcore.size());

  std::vector<std::pair<unsigned, Term>> score_vec;
  for (const auto & t : unsatcore) {
    auto score = GetScore(t);
    score_vec.push_back(std::make_pair(score, t));
  }
  std::sort(score_vec.begin(), score_vec.end());
  std::unordered_map<Term, unsigned> term_id_map;
  TermList sorted_unsatcore;
  unsigned pidx = 0;
  for (auto pos = score_vec.rbegin(); pos != score_vec.rend(); ++pos) {
    sorted_unsatcore.push_back(pos->second);
    D(4, "{}: Score: {} expr: {}", pidx, pos->first, pos->second->to_string());
    term_id_map.emplace(pos->second, pidx++);
  }
  push_solver_context();
  disable_all_labels();
  syntax_analysis::reduce_unsat_core_linear(base, sorted_unsatcore, solver_);
  pop_solver_context();

#ifdef DEBUG
  std::cout << "Total : " << preds_nxt.size() << " {";
  for (const auto & p : sorted_unsatcore) {
    auto pos = term_id_map.find(p);
    assert(pos != term_id_map.end());
    std::cout << pos->second << ",";
  }
  std::cout << "}" << std::endl;
#endif
  assert(!sorted_unsatcore.empty());

  TermVec curr_lits;
  curr_lits.reserve(sorted_unsatcore.size());
  for (const auto & l : sorted_unsatcore) {
    curr_lits.push_back(ts_.curr(l));  // make sure we also convert input back
  }
  return ic3formula_negate(ic3formula_conjunction(curr_lits));
}  // select_predicates_btor

IC3Formula SygusPdr::select_predicates_generic(const Term & base,
                                               const TermVec & preds_nxt)
{
  TermVec unsatcore;
  TermVec sorted_unsatcore;
  D(3, "[IterativeReduction] Initial core size {}", preds_nxt.size());

  UnorderedTermSet pred_set(preds_nxt.begin(), preds_nxt.end());
  TermVec deduplicate_preds_nxt(pred_set.begin(), pred_set.end());

  reducer_.reduce_assump_unsatcore(
      base, deduplicate_preds_nxt, unsatcore, NULL, 0, options_.random_seed_);

  D(3, "[IterativeReduction] End core size {}", unsatcore.size());
  std::vector<std::pair<unsigned, Term>> score_vec;
  for (const auto & t : unsatcore) {
    auto score = GetScore(t);
    score_vec.push_back(std::make_pair(score, t));
  }
  std::sort(score_vec.begin(), score_vec.end());
  std::unordered_map<Term, unsigned> term_id_map;
  unsigned pidx = 0;
  for (auto pos = score_vec.rbegin(); pos != score_vec.rend(); ++pos) {
    sorted_unsatcore.push_back(pos->second);
    D(4, "{}: Score: {} expr: {}", pidx, pos->first, pos->second->to_string());
    term_id_map.emplace(pos->second, pidx++);
  }
  unsatcore.clear();  // reuse this
  reducer_.linear_reduce_assump_unsatcore(
      base, sorted_unsatcore, unsatcore, NULL, 0);
#ifdef DEBUG
  std::cout << "Total : " << preds_nxt.size() << " {";
  for (const auto & p : unsatcore) {
    auto pos = term_id_map.find(p);
    assert(pos != term_id_map.end());
    std::cout << pos->second << ",";
  }
  std::cout << "}" << std::endl;
#endif
  assert(!unsatcore.empty());

  TermVec curr_lits;
  curr_lits.reserve(unsatcore.size());
  for (const auto & l : unsatcore) {
    curr_lits.push_back(ts_.curr(l));  // make sure we also convert input back
  }
  return ic3formula_negate(ic3formula_conjunction(curr_lits));
}  // select_predicates

bool SygusPdr::propose_new_terms(syntax_analysis::IC3FormulaModel * pre_model,
                                 syntax_analysis::IC3FormulaModel * post_model,
                                 const smt::Term & F_T_not_cex,
                                 const smt::Term & Init_prime,
                                 bool failed_at_init)
{
  unsigned proposing_new_terms_round = 0;
  unsigned n_new_terms;

  do {
    n_new_terms = sygus_term_manager_.GetMoreTerms(pre_model,
                                                   post_model,
                                                   *(term_learner_.get()),
                                                   ts_.trans(),
                                                   failed_at_init,
                                                   options_.sygus_term_mode_);
    D(3,
      "[propose-new-term] Round {}. Get {} new terms.",
      proposing_new_terms_round,
      n_new_terms);
    if (n_new_terms != 0) {
      syntax_analysis::PredConstructor pred_collector_(
          to_next_func_,
          solver_,
          post_model,  // IC3FormulaModel * cex
          setup_cex_info(post_model)
          // you need to set this up every time you use it
          // this is needed to eval new terms
      );
      Term base = solver_->make_term(
          Or,
          solver_->make_term(And, F_T_not_cex, pre_model->to_expr()),
          Init_prime);
      push_solver_context();
      disable_all_labels();
      solver_->assert_formula(base);
      const auto & pred_nxt = pred_collector_.GetAllPredNext();
      for (const auto & p : pred_nxt) solver_->assert_formula(p);
      auto res = check_sat();
      pop_solver_context();

      if (res.is_unsat()) return true;
    }
  } while (n_new_terms != 0);
  // at this point no new terms
  return false;
}  // propose_new_terms

syntax_analysis::PerCexInfo & SygusPdr::setup_cex_info(
    syntax_analysis::IC3FormulaModel * post_model)
{
  assert(post_model);
  auto cex_term_map_pos = cex_term_map_.find(post_model);
  if (cex_term_map_pos == cex_term_map_.end()) {
    bool nouse;
    std::tie(cex_term_map_pos, nouse) = cex_term_map_.emplace(
        post_model,
        syntax_analysis::PerCexInfo(
            sygus_term_manager_.GetAllTermsForVarsInModel(
                post_model,
                solver_,
                options_.sygus_term_mode_,
                options_.sygus_term_extract_depth_,
                options_.sygus_initial_term_width_,
                options_.sygus_initial_term_inc_,
                options_.sygus_accumulated_term_bound_),
            op_uf_assumptions_.size() /*constraint_count*/));
  }  // if no terms, set up them first

  syntax_analysis::PerCexInfo & per_cex_info = cex_term_map_pos->second;

  // then evaluate the terms on the cex
  push_solver_context();
  disable_all_labels();
  for (const auto & c : op_uf_assumptions_) solver_->assert_formula(c);
  solver_->assert_formula(post_model->to_expr());
  auto res = check_sat();
  assert(res.is_sat());

  bool reset_due_to_more_refinement =
      per_cex_info.prev_refine_constraint_count < op_uf_assumptions_.size();
  per_cex_info.prev_refine_constraint_count = op_uf_assumptions_.size();

  // for each width
  for (const auto & width_term_const_pair : per_cex_info.varset_info.terms) {
    auto width = width_term_const_pair.first;
    if (reset_due_to_more_refinement)
      per_cex_info.prev_per_width_term_num[width].term_num = 0;
    // when to update this ? after predicates are generated
    unsigned nc = per_cex_info.prev_per_width_term_num[width].const_num;
    unsigned nt = per_cex_info.prev_per_width_term_num[width].term_num;

    auto nt_end = width_term_const_pair.second.terms.size();
    auto nc_end = width_term_const_pair.second.constants.size();
    // cache the terms and constants value under the cex
    for (unsigned tidx = nt; tidx < nt_end; ++tidx) {
      const auto & t = width_term_const_pair.second.terms.at(tidx);
      per_cex_info.terms_val_under_cex.emplace(
          t, syntax_analysis::eval_val(solver_->get_value(t)->to_string()));
    }
    for (unsigned cidx = nc; cidx < nc_end; ++cidx) {
      const auto & c = width_term_const_pair.second.constants.at(cidx);
      per_cex_info.terms_val_under_cex.emplace(c, c->to_string());
    }
  }  // eval terms on cex
  pop_solver_context();

  if (reset_due_to_more_refinement) per_cex_info.ResetPredicates();

  return per_cex_info;
}  // setup_cex_info

std::pair<IC3Formula, syntax_analysis::IC3FormulaModel *>
SygusPdr::ExtractInitPrimeModel(const Term & p_prime)
{
  UnorderedTermSet varlist;
  partial_model_getter.GetVarList(p_prime, varlist);

  Term conj_partial;
  TermVec conjvec_partial;
  syntax_analysis::IC3FormulaModel::cube_t cube_partial;

  for (const auto & v : varlist) {
    Term val = solver_->get_value(v);
    Term curr_v = ts_.curr(v);
    auto eq = solver_->make_term(Op(PrimOp::Equal), curr_v, val);
    cube_partial.emplace(curr_v, val);
    conjvec_partial.push_back(eq);
    if (conj_partial) {
      conj_partial = solver_->make_term(Op(PrimOp::And), conj_partial, eq);
    } else {
      conj_partial = eq;
    }
  }
  if (conj_partial == nullptr) {
    conj_partial = solver_true_;
    assert(conjvec_partial.empty());
    conjvec_partial.push_back(solver_true_);
  }

  syntax_analysis::IC3FormulaModel * partial_model =
      new syntax_analysis::IC3FormulaModel(std::move(cube_partial),
                                           conj_partial);

  assert(partial_model);
  model2cube_.emplace(conj_partial, partial_model);

  return std::make_pair(
      IC3Formula(conj_partial, { conj_partial }, false /*not a disjunction*/),
      partial_model);
}  // ExtractInitPrimeModel

std::tuple<IC3Formula,
           syntax_analysis::IC3FormulaModel *,
           syntax_analysis::IC3FormulaModel *>
SygusPdr::ExtractPartialAndFullModel(
    syntax_analysis::IC3FormulaModel * post_model)
{
  // we get varset from post_model, this will save us from one term walk
  UnorderedTermSet varset;
  post_model->get_varset(varset);
  TermVec pre_var;
  for (const auto & v : varset) {
    if (IN(v, ts_.statevars()) && !IN(v, no_next_vars_)) {
      auto update_pos = ts_.state_updates().find(v);
      assert(update_pos != ts_.state_updates().end());
      pre_var.push_back(update_pos->second);
    }  // if it is input then we should not need
  }
  UnorderedTermSet varlist;
  partial_model_getter.GetVarListForAsts(pre_var, varlist);

  // get varset, convert to current var (take care of input vars)
  // extract var from it, with/without inputs
  // and also extract the full model
  {
    Term conj_partial;
    TermVec conjvec_partial;
    syntax_analysis::IC3FormulaModel::cube_t cube_partial;

    Term conj_full;
    syntax_analysis::IC3FormulaModel::cube_t cube_full;

    for (const auto & v : varlist) {
      Term val = solver_->get_value(v);
      auto eq = solver_->make_term(Op(PrimOp::Equal), v, val);

      cube_full.emplace(v, val);
      if (conj_full)
        conj_full = solver_->make_term(Op(PrimOp::And), conj_full, eq);
      else
        conj_full = eq;

      if (keep_var_in_partial_model(v)) {
        cube_partial.emplace(v, val);
        conjvec_partial.push_back(eq);
        if (conj_partial) {
          conj_partial = solver_->make_term(Op(PrimOp::And), conj_partial, eq);
        } else {
          conj_partial = eq;
        }
      }  // end of partial model
    }
    assert(conj_full);  // at least there should be var in the model
    if (conj_full == nullptr) conj_full = solver_true_;

    if (conj_partial == nullptr) {
      conj_partial = solver_true_;
      assert(conjvec_partial.empty());
      conjvec_partial.push_back(solver_true_);
    }
    // this is called from inductive_gen which should not create must block goal
    syntax_analysis::IC3FormulaModel * partial_model =
        new syntax_analysis::IC3FormulaModel(std::move(cube_partial),
                                             conj_partial);

    syntax_analysis::IC3FormulaModel * full_model =
        new syntax_analysis::IC3FormulaModel(std::move(cube_full), conj_full);

    assert(partial_model && full_model);
    model2cube_.emplace(conj_partial, partial_model);

    return std::make_tuple(
        IC3Formula(conj_partial, { conj_partial }, false /*not a disjunction*/),
        partial_model,
        full_model);
  }  // end model making block
}  // ExtractPartialAndFullModel

std::pair<IC3Formula, syntax_analysis::IC3FormulaModel *>
SygusPdr::ExtractPartialModel(const Term & p)
{
  // extract using keep_var_in_partial_model
  assert(ts_.no_next(p));

  UnorderedTermSet varlist;
  Term bad_state_no_nxt = next_curr_replace(ts_.next(p));

  // we need to make sure input vars are mapped to next input vars

  D(4, "[PartialModel] prime state : {}", bad_state_no_nxt->to_string());
  if (has_assumptions) {
    D(4,
      "[PartialModel] assumptions (mapped): {}",
      constraints_curr_var_.size());
    unsigned idx = 0;
    for (const auto & c : constraints_curr_var_)
      D(4, "[PartialModel] assumption #{} : {}", idx++, c->to_string());
    constraints_curr_var_.push_back(bad_state_no_nxt);
    partial_model_getter.GetVarListForAsts(constraints_curr_var_, varlist);
    constraints_curr_var_.pop_back();
  } else {
    partial_model_getter.GetVarList(bad_state_no_nxt, varlist);
  }

  {
    D(4, "[PartialModel] before cutting vars: ");
    for (const auto & v : varlist)
      D(4,
        "[PartialModel] {} := {} ",
        v->to_string(),
        solver_->get_value(v)->to_string());
    D(4, "[PartialModel] ------------------- ");
  }

  Term conj_partial;
  TermVec conjvec_partial;
  syntax_analysis::IC3FormulaModel::cube_t cube_partial;

  for (const auto & v : varlist) {
    Term val = solver_->get_value(v);
    auto eq = solver_->make_term(Op(PrimOp::Equal), v, val);
    if (keep_var_in_partial_model(v)) {
      cube_partial.emplace(v, val);
      conjvec_partial.push_back(eq);
      if (conj_partial) {
        conj_partial = solver_->make_term(Op(PrimOp::And), conj_partial, eq);
      } else {
        conj_partial = eq;
      }
    }  // end of partial model
  }
  if (conj_partial == nullptr) {
    conj_partial = solver_true_;
    assert(conjvec_partial.empty());
    conjvec_partial.push_back(solver_true_);
  }
  syntax_analysis::IC3FormulaModel * partial_model =
      new syntax_analysis::IC3FormulaModel(std::move(cube_partial),
                                           conj_partial);

  assert(partial_model);
  model2cube_.emplace(conj_partial, partial_model);

  return std::make_pair(
      IC3Formula(conj_partial, { conj_partial }, false /*not a disjunction*/),
      partial_model);
}  // SygusPdr::ExtractPartialAndFullModel

void SygusPdr::predecessor_generalization(size_t i,
                                          const Term & cterm,
                                          IC3Formula & pred)
{
  // used in rel_ind_check
  // extract the model based on var list
  // NOTE: i may be incorrect, it is given as F/\T->(here is i)

  // no need to pop (pop in rel_ind_check)
  // return the model and build IC3FormulaModel
  auto partial_full_model = ExtractPartialModel(cterm);
  pred = partial_full_model.first;
}  // generalize_predecessor

bool SygusPdr::keep_var_in_partial_model(const Term & v) const
{
  if (options_.sygus_use_operator_abstraction_) {
    assert(op_abstractor_ != nullptr);
    const auto & dummy_inputs = op_abstractor_->dummy_inputs();
    if (dummy_inputs.find(v) != dummy_inputs.end()) return true;
  }

  if (has_assumptions) {  // must keep input vars
    return (ts_.is_curr_var(v));
  }

  return ts_.is_curr_var(v) && !IN(v, no_next_vars_);
}  // keep_var_in_partial_model

void SygusPdr::abstract()
{
  // called in initialize()
}  // SygusPdr::abstract()

RefineResult SygusPdr::refine()
{
  if (!options_.sygus_use_operator_abstraction_) return REFINE_NONE;

  assert(witness_length() > 0);

  TermVec new_constraints;
  bool succ =
      op_abstractor_->refine_with_constraints(cex_, bad_, new_constraints);
  assert(succ == (!new_constraints.empty()));
  if (succ) {
    for (const auto & c : new_constraints) {
      op_uf_assumptions_.push_back(c);
      ts_.add_constraint(c, true);
      D(4, "[Refine] : {} ", c->to_string());

      // rework the constraints buffers
      assert(ts_.no_next(c));
      smt::UnorderedTermSet vars;
      get_free_symbolic_consts(c, vars);
      if (!vars.empty()) {
        constraints_curr_var_.push_back(c);
        // translate input_var to next input_var
        // but the state var ...
        constraints_curr_var_.push_back(next_curr_replace(ts_.next(c)));
        has_assumptions = true;  // previously it can also be true
      }
    }
    reset_solver();  // because tr/init could be different now

    return REFINE_SUCCESS;
  }  // else
  return REFINE_NONE;
}  // SygusPdr::refine()

// ------------------------------------------------------------------
//  May Block

bool SygusPdr::try_recursive_block_goal(const IC3Formula & to_block,
                                        unsigned fidx)
{
  assert(!solver_context_);
  ProofGoalQueue proof_goals;
  proof_goals.new_proof_goal(to_block, fidx, nullptr);
  while (!proof_goals.empty()) {
    const ProofGoal * pg = proof_goals.top();
    if (pg->idx == 0) {  // fail
      // do we refine ? no, because this is just a may proof goal
      return false;
    }

    if (is_blocked(pg)) {
      assert(pg == proof_goals.top());
      proof_goals.pop();
      continue;
    }

    //  try and see if it is blockable
    IC3Formula
        collateral;  // populated by rel_ind_check: if unsat collateral:=target
                     // if sat collateral:=partial model
    if (rel_ind_check_may_block(pg->idx, pg->target, collateral)) {
      // this proof goal can be blocked
      assert(!solver_context_);
      assert(collateral.term);
      logger.log(3, "Blocking term at frame {}: {}", pg->idx, pg->target.term);

      // remove the proof goal now that it has been blocked
      assert(pg == proof_goals.top());
      proof_goals.pop();
      assert(collateral.term == pg->target.term);
      collateral = inductive_generalization(pg->idx, collateral);

      size_t idx = find_highest_frame(pg->idx, collateral);
      assert(idx >= pg->idx);

      assert(collateral.disjunction);
      assert(collateral.term);
      assert(collateral.children.size());
      SygusPdr::constrain_frame(idx, collateral);
      logger.log(3,
                 "Blocking term at frame {}: {} , \n --- with {}",
                 pg->idx,
                 pg->target.term,
                 collateral.term->to_string());
    } else {
      if (collateral.term == solver_true_)
        return false;  // when we intersect with 0
      // could not block this proof goal
      assert(collateral.term);
      proof_goals.new_proof_goal(collateral, pg->idx - 1, pg);
    }
  }  // end of while
  assert(proof_goals.empty());
  return true;
}  // try_recursive_block_goal_at_or_before

// differences from the ic3base class
// (a) no -c (b) when unsat, out = c (c) new check
bool SygusPdr::rel_ind_check_may_block(size_t i,
                                       const IC3Formula & c,
                                       IC3Formula & out)
{
  assert(i > 0);
  assert(i < frames_.size());
  // expecting to be the polarity for proof goals, not frames
  // e.g. a conjunction
  assert(!c.disjunction);

  assert(solver_context_ == 0);
  push_solver_context();

  // F[i-1]
  assert_frame_labels(i - 1);
  // -c
  // solver_->assert_formula(solver_->make_term(Not, c.term));
  // Trans
  assert_trans_label();

  solver_->assert_formula(ts_.next(c.term));

  Result r = check_sat();
  if (r.is_sat()) {
    // get model
    predecessor_generalization(i, c.term, out);
    assert(out.term);
    assert(out.children.size());
    assert(!out.disjunction);  // expecting a conjunction
    assert(ic3formula_check_valid(out));

    pop_solver_context();
    return r.is_unsat();
  }
  pop_solver_context();
  push_solver_context();
  solver_->assert_formula(init_label_);
  solver_->assert_formula(c.term);
  r = check_sat();
  if (r.is_sat()) {
    // get model at 0
    out.term = solver_true_;
    out.children = { solver_true_ };
    out.disjunction = false;
    pop_solver_context();
    return r.is_unsat();
  }

  out = c;
  assert(!out.disjunction);  // we need it to be a model
  pop_solver_context();
  assert(!r.is_unknown());
  return r.is_unsat();
}  // rel_ind_check_may_block

// ---------------------------------------------------------------------
//   Override for: bad = model(F/\T -> P')
// ---------------------------------------------------------------------

void SygusPdr::disable_all_labels()
{
#define NOT(x) (solver_->make_term(Not, (x)))
  solver_->assert_formula(NOT(init_label_));
  solver_->assert_formula(NOT(trans_label_));
  solver_->assert_formula(NOT(bad_label_));
  for (const auto & fl : frame_labels_) solver_->assert_formula(NOT(fl));
#undef NOT
}

}  // namespace pono
