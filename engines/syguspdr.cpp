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
#include "utils/logger.h"
#include "utils/sygus_predicate_constructor.h"

using namespace smt;

namespace pono {

void CustomFunctionalTransitionSystem::make_nextvar_for_inputs() {
  // make next variables for the input variables
  for (const auto & i : inputvars_) {
    const auto & input_name = term_to_name_.at(i);
    const auto next_input_name = input_name + ".next";
    const auto sort = i->get_sort();
    auto next_inputvar = solver_->make_symbol(next_input_name, sort);

    next_inputvars_.emplace(i,next_inputvar);
    next_map_[i] = next_inputvar;
    curr_map_[next_inputvar] = i;
    name_term(next_inputvar->to_string(), next_inputvar);
  }
}

// ----------------------------------------------------------------
// Helper functions
// ----------------------------------------------------------------

unsigned SygusPdr::GetScore(const Term & t) {
  const auto & score_map = term_score_walker_.GetScoreMap();
  auto pos = score_map.find(t);
  if (pos != score_map.end())
    return pos->second.score;
  term_score_walker_.WalkBFS(t);
  pos = score_map.find(t);
  assert (pos != score_map.end());
  return pos->second.score;
}


// ----------------------------------------------------------------
// Constructor & Destructor
// ----------------------------------------------------------------
SygusPdr::SygusPdr(Property & p, const SmtSolver & slv,
                             PonoOptions opt)
  : super(p, slv, opt),
    partial_model_getter(solver_),
    custom_ts_(NULL),
    has_assumptions(true) // most conservative way
{
  solver_->set_opt("produce-unsat-cores", "true");
}

// destruct -- free the memory
SygusPdr::~SygusPdr() {
  for (const auto & ic3f_mptr : model2cube_) {
    if(ic3f_mptr.second)
      delete ic3f_mptr.second;
  }
  
  //  for (const auto & partial2full : to_full_model_map_) {
  //    if (partial2full.second)
  //      delete partial2full.second; // delete the full models
  //  }

  if (custom_ts_)
    delete custom_ts_;
} // SygusPdr::~SygusPdr

// ----------------------------------------------------------------
// IC3 override
// ----------------------------------------------------------------


// init 
void SygusPdr::initialize()
{
  
  if (initialized_)
    return;

  custom_ts_ = 
    new CustomFunctionalTransitionSystem(orig_ts_, TermTranslator( orig_ts_.solver()) );
  ts_ = custom_ts_;
  custom_ts_->make_nextvar_for_inputs();

  bad_next_ = custom_ts_->next(bad_);

  // has_assumption
  has_assumptions = ! (custom_ts_->constraints().empty());
  for (const auto & c : custom_ts_->constraints()) {
    if (custom_ts_->no_next(c)) {
      constraints_curr_var_.push_back(c);
      // translate input_var to next input_var
      // but the state var ...
      constraints_curr_var_.push_back(
        custom_ts_->to_next_func(
          solver_->substitute(c,
            custom_ts_->input_var_to_next_map())));
    } // else skip
  }

  super::initialize();
  { // add P to F[1]
    auto prop = smart_not(bad_);
    constrain_frame(1, IC3Formula(prop, {prop}, false), true);
  }
  
  // initialize the caches  
  // extract the operators
  op_extract_ = std::make_unique<syntax_analysis::OpExtractor>();
  op_extract_->WalkBFS(custom_ts_->init());
  op_extract_->WalkBFS(custom_ts_->trans());
  op_extract_->GetSyntaxConstruct().RemoveConcat();
  op_extract_->GetSyntaxConstruct().RemoveExtract();
  op_extract_->GetSyntaxConstruct().AndOrConvert();
  op_extract_->GetSyntaxConstruct().RemoveUnusedStructure();
    

  { // 1. register terms to find exprs
    // 2. extract parent from the same terms
    for (auto && v_nxtexpr_pair : custom_ts_->state_updates()) {
      assert(custom_ts_->no_next(v_nxtexpr_pair.second));
      sygus_term_manager_.RegisterTermsToWalk(v_nxtexpr_pair.second);
      parent_of_terms_.WalkBFS(v_nxtexpr_pair.second);
    }
    sygus_term_manager_.RegisterTermsToWalk(custom_ts_->init());
    parent_of_terms_.WalkBFS(custom_ts_->init());

    for (const auto & c : custom_ts_->constraints()) {
      if (custom_ts_->no_next(c)) {
        sygus_term_manager_.RegisterTermsToWalk(c);
        parent_of_terms_.WalkBFS(c);
      }
    }

    sygus_term_manager_.RegisterTermsToWalk(property_.prop());
    parent_of_terms_.WalkBFS(property_.prop());
  }

  // cache two lambda functions for sygus enum
  to_next_func_ = [&] (const Term & v) -> Term {
    return custom_ts_->next(v);
  };

  { // now create TermLearner
    term_learner_.reset(new syntax_analysis::TermLearner(
      to_next_func_, solver_,  // you need to make trans dynamic
      parent_of_terms_
    ));
  }

} // initialize


void SygusPdr::check_ts() const
{
  // check if there are arrays or uninterpreted sorts and fail if so
  if (!ts_->is_functional())
    throw PonoException(
      "SyGuS PDR only supports functional transition systems.");
    // check if there are arrays or uninterpreted sorts and fail if so
  for (auto vec : { ts_->statevars(), ts_->inputvars() }) {
    for (auto st : vec) {
      SortKind sk = st->get_sort()->get_sort_kind();
      if (sk == ARRAY) {
        throw PonoException("SyGuS PDR does not support arrays yet");
      } else if (sk == UNINTERPRETED) {
        throw PonoException(
            "SyGuS PDR does not support uninterpreted sorts yet.");
      }
    }
  }
} // check_ts




bool SygusPdr::ic3formula_check_valid(const IC3Formula & u) const
{
  const Sort &boolsort = solver_->make_sort(BOOL);
  // check that children are literals
  Term pred;
  Op op;
  for (const auto &c : u.children) {
    if (c->get_sort() != boolsort) {
      logger.log(3, "ERROR SyGuS PDR IC3Formula contains non-boolean atom: {}", c);
      return false;
    }
  }
  // got through all checks without failing
  return true;
}


IC3Formula SygusPdr::get_model_ic3formula(
    TermVec * out_inputs,
    TermVec * out_nexts) const {
  assert(false); // should not use this
  throw PonoException(
            "SyGuS PDR must be used with generalize_predecessor.");
  return IC3Formula();
} // SygusPdr::get_model_ic3formula


std::vector<IC3Formula> SygusPdr::inductive_generalization(
    size_t i, const IC3Formula & c) {
  assert(i>0);
  // this is to produce the lemma
  // find the model
  auto model_pos = model2cube_.find( c.term );
  assert(model_pos != model2cube_.end());
  syntax_analysis::IC3FormulaModel * post_model = model_pos->second;
  syntax_analysis::PerCexInfo & per_cex_info = setup_cex_info(post_model);
  // find the related info
  
  syntax_analysis::PredConstructor pred_collector_(
    to_next_func_, solver_,
    post_model, // IC3FormulaModel * cex 
    per_cex_info
  );

  Term Fprev = get_frame_term(i - 1);
  Term T = ts_->trans();
  Term Init_prime = ts_->next(ts_->init());
  const Term & cex = c.term;
  Term not_cex = smart_not(cex);
  Term cex_prime = ts_->next(cex);
  Term F_T_not_cex = make_and( {Fprev, not_cex, T} );
  Term base = 
    solver_->make_term(Or,
      F_T_not_cex,
      Init_prime);


  IC3Formula pre_formula;
  syntax_analysis::IC3FormulaModel * pre_model = NULL;
  bool failed_at_init;

  bool insufficient_pred = false;
  do {
    { // step 1 - check if the pred are good if not construct the model
      const auto & pred_nxt_ = pred_collector_.GetAllPredNext();
      push_solver_context();
      disable_all_labels();
      // check ( F[i-1] /\ NOT(c) /\ T ) \/ init' AND pred_next sat

      solver_->assert_formula(base);
      for (const auto & pred_ : pred_nxt_)
        solver_->assert_formula(pred_);
      
      Result r = check_sat();
      if (r.is_sat()) {
        // TODO: extract model and do recursive block
        // pre_model will be put into the map and free at class destructor

        insufficient_pred = true;
        failed_at_init = solver_->get_value(Init_prime)->to_int() != 0;
        if (!failed_at_init) {
          std::tie(pre_formula, pre_model) =
            ExtractPartialModel( cex );
        } else {
          std::tie(pre_formula, pre_model) =
            ExtractInitPrimeModel( Init_prime ); // this can be only the prime vars
        }
        // note Here you must give a formula with current variables
      } else
        insufficient_pred = false;
      pop_solver_context();
    }

    if (insufficient_pred) {
      // extract Proof goal from partial_model
      if (!failed_at_init) {
        if(try_recursive_block_goal_at_or_before(pre_formula, i-1))
          continue; // see if we have next that we may need to block
      } // if we failed at init, then we will anyway need to 
      // propose new terms

      bool succ = 
        propose_new_terms(pre_model, post_model, per_cex_info,
          F_T_not_cex, Init_prime, failed_at_init);
      assert(succ);
      // it has loop inside
      // and it must succeed because eventually we will end with bit-level things
    }
  }while(insufficient_pred);

  // at this point we have enough preds
  // reduce the preds
  return {select_predicates(
      base,
      pred_collector_.GetAllPredNext()
    )};
} // inductive_generalization

IC3Formula SygusPdr::select_predicates(const Term & base, const TermVec & preds_nxt) {
  TermVec unsatcore;
  TermVec sorted_unsatcore;
  reducer_.reduce_assump_unsatcore(base, preds_nxt, unsatcore, NULL, 0, options_.random_seed_);
  std::vector<std::pair<unsigned, Term>> score_vec;
  for (const auto & t : unsatcore) {
    auto score = GetScore(t);
    score_vec.push_back(std::make_pair(score, t));
  }
  std::sort(score_vec.begin(),score_vec.end());
  for (auto pos = score_vec.rbegin(); pos != score_vec.rend(); ++ pos) {
    sorted_unsatcore.push_back(pos->second);
  }
  unsatcore.clear(); // reuse this
  reducer_.linear_reduce_assump_unsatcore(base, sorted_unsatcore, unsatcore, NULL, 0);
  
  TermVec curr_lits;
  curr_lits.reserve(unsatcore.size());
  for (const auto &l : unsatcore) 
    curr_lits.push_back(ts_->curr(l));
  return ic3formula_negate(ic3formula_conjunction(curr_lits));
} // select_predicates

bool SygusPdr::propose_new_terms(
    syntax_analysis::IC3FormulaModel * pre_model,
    syntax_analysis::IC3FormulaModel * post_model,
    syntax_analysis::PerCexInfo & per_cex_info,
    const smt::Term & F_T_not_cex,
    const smt::Term & Init_prime,
    bool failed_at_init
    ) {
  
  unsigned proposing_new_terms_round = 0;
  unsigned n_new_terms;

  do {
    n_new_terms = 
      sygus_term_manager_.GetMoreTerms(
        pre_model, post_model, *(term_learner_.get()),
        ts_->trans(),
        failed_at_init);
    logger.log(3, "[propose-new-term] Round {}. Get {} new terms.", proposing_new_terms_round, n_new_terms);
    if (n_new_terms != 0) {
      syntax_analysis::PredConstructor pred_collector_(
        to_next_func_, solver_,
        post_model, // IC3FormulaModel * cex 
        per_cex_info
      );
      Term base = 
        solver_->make_term(Or,
          solver_->make_term(And, F_T_not_cex, pre_model->to_expr()),
          Init_prime);
      push_solver_context();
      disable_all_labels();
        solver_->assert_formula(base);
        const auto & pred_nxt = pred_collector_.GetAllPredNext();
        for (const auto & p: pred_nxt)
          solver_->assert_formula(p);
      auto res = solver_->check_sat();
      pop_solver_context();

      if (res.is_unsat())
        return true;
    }
  } while(n_new_terms != 0);
  // at this point no new terms
  return false;
} // propose_new_terms

syntax_analysis::PerCexInfo & SygusPdr::setup_cex_info (syntax_analysis::IC3FormulaModel * post_model) {
  assert(post_model);
  auto cex_term_map_pos = cex_term_map_.find(post_model);
  if (cex_term_map_pos == cex_term_map_.end()) {
    bool nouse;
    std::tie(cex_term_map_pos ,  nouse) =
      cex_term_map_.emplace(post_model, 
          syntax_analysis::PerCexInfo( 
            sygus_term_manager_.GetAllTermsForVarsInModel(post_model, solver_, options_.term_option ) ) );
  } // if no terms, set up them first

  syntax_analysis::PerCexInfo & per_cex_info = cex_term_map_pos->second;
  unsigned nterm = 0;

  // then evalute the terms on the cex
  push_solver_context();
  disable_all_labels();

  solver_->assert_formula( post_model->to_expr() );
  auto res = solver_->check_sat();
  assert (res.is_sat());
  // for each witdh
  for (const auto & width_term_const_pair : per_cex_info.varset_info.terms) {
    auto width = width_term_const_pair.first;
    unsigned nc = per_cex_info.prev_per_width_term_num[width].const_num; // when to update this
    unsigned nt = per_cex_info.prev_per_width_term_num[width].term_num;
    auto nt_end = width_term_const_pair.second.terms.size();
    auto nc_end = width_term_const_pair.second.constants.size();
    // cache the terms and constants value under the cex
    for (unsigned tidx = nt; tidx < nt_end ; ++tidx) {
      const auto & t = width_term_const_pair.second.terms.at(tidx);
      per_cex_info.terms_val_under_cex.emplace(
        t, syntax_analysis::eval_val( solver_->get_value(t)->to_string() ));
      ++ nterm;
    }
    for (unsigned cidx = nc ; cidx <  nc_end; ++cidx) {
      const auto & c = width_term_const_pair.second.constants.at(cidx);
      per_cex_info.terms_val_under_cex.emplace(
        c, c->to_string() );
      ++ nterm;
    }
  } // eval terms on cex
  pop_solver_context();
  return per_cex_info;
} // setup_cex_info

IC3Formula SygusPdr::get_initial_bad_model() {
  Term conj_full;
  for (const auto & v : ts_->statevars()) {
    Term val = solver_->get_value(v);
    auto eq = solver_->make_term(Op(PrimOp::Equal), v,val );
    
    if (conj_full)
      conj_full = solver_->make_term(Op(PrimOp::And), conj_full, eq);
    else
      conj_full = eq;
  }

  for (const auto & v : ts_->inputvars()) {
    Term val = solver_->get_value(v);
    auto eq = solver_->make_term(Op(PrimOp::Equal), v,val );
    
    if (conj_full)
      conj_full = solver_->make_term(Op(PrimOp::And), conj_full, eq);
    else
      conj_full = eq;
  }

  return IC3Formula(conj_full, {conj_full}, false /*not a disjunction*/ );
} // get_initial_bad_model

std::pair<IC3Formula, syntax_analysis::IC3FormulaModel *>
  SygusPdr::ExtractInitPrimeModel(const Term & p_prime) {

  UnorderedTermSet varlist;
  partial_model_getter.GetVarList(p_prime, varlist);

  Term conj_partial;
  TermVec conjvec_partial;
  syntax_analysis::IC3FormulaModel::cube_t cube_partial;

  
  for (const auto & v : varlist) {
    Term val = solver_->get_value(v);
    auto eq = solver_->make_term(Op(PrimOp::Equal), v,val );
    if (keep_var_in_partial_model(v)) {
      cube_partial.emplace(v,val);
      conjvec_partial.push_back( eq );
      if (conj_partial) {
        conj_partial = solver_->make_term(Op(PrimOp::And), conj_partial, eq);
      } else {
        conj_partial = eq;
      }
    } // end of partial model
  }
  syntax_analysis::IC3FormulaModel * partial_model = 
    new syntax_analysis::IC3FormulaModel(std::move(cube_partial), conj_partial);

  assert(partial_model);
  model2cube_.emplace(conj_partial, partial_model);

  return std::make_pair(
    IC3Formula(conj_partial, conjvec_partial, false /*not a disjunction*/ ),
    partial_model
  );
} // ExtractInitPrimeModel 

std::pair<IC3Formula, syntax_analysis::IC3FormulaModel *>
    SygusPdr::ExtractPartialModel(const Term & p) {
  // extract using keep_var_in_partial_model  
  assert(custom_ts_->no_next(p));

  UnorderedTermSet varlist;
  Term bad_state_no_nxt = custom_ts_->to_next_func(p);
  if (has_assumptions) {
    constraints_curr_var_.push_back(bad_state_no_nxt);
    partial_model_getter.GetVarListForAsts(constraints_curr_var_, varlist);
    constraints_curr_var_.pop_back();
  } else {
    partial_model_getter.GetVarList(bad_state_no_nxt, varlist);
  }

  Term conj_partial;
  TermVec conjvec_partial;
  syntax_analysis::IC3FormulaModel::cube_t cube_partial;

  
  for (const auto & v : varlist) {
    Term val = solver_->get_value(v);
    auto eq = solver_->make_term(Op(PrimOp::Equal), v,val );
    if (keep_var_in_partial_model(v)) {
      cube_partial.emplace(v,val);
      conjvec_partial.push_back( eq );
      if (conj_partial) {
        conj_partial = solver_->make_term(Op(PrimOp::And), conj_partial, eq);
      } else {
        conj_partial = eq;
      }
    } // end of partial model
  }
  syntax_analysis::IC3FormulaModel * partial_model = 
    new syntax_analysis::IC3FormulaModel(std::move(cube_partial), conj_partial);

  assert(partial_model);
  model2cube_.emplace(conj_partial, partial_model);

  return std::make_pair(
    IC3Formula(conj_partial, conjvec_partial, false /*not a disjunction*/ ),
    partial_model
  );
} // SygusPdr::ExtractPartialAndFullModel


IC3Formula SygusPdr::generalize_predecessor(size_t i, const IC3Formula & c) {
  // used in rel_ind_check
  // extract the model based on var list

  // no need to pop (pop in rel_ind_check)
  // return the model and build IC3FormulaModel
  auto partial_full_model = ExtractPartialModel(c.term);
  return partial_full_model.first;
} // generalize_predecessor


bool SygusPdr::keep_var_in_partial_model(const Term & v) const {
  if (has_assumptions) { // must keep input vars
    if(custom_ts_->is_curr_var(v) || custom_ts_->inputvars().find(v) != custom_ts_->inputvars().end() )
      return true;
    return false;
  }
  return custom_ts_->is_curr_var(v);
} // keep_var_in_partial_model



void SygusPdr::abstract() {
  // called in initialize()
} // SygusPdr::abstract()


RefineResult SygusPdr::refine() {
  return REFINE_NONE;
} // SygusPdr::refine()


bool SygusPdr::try_recursive_block_goal_at_or_before(const IC3Formula & to_block, unsigned fidx) {
  add_proof_goal(to_block, fidx, NULL);

  while (has_proof_goals()) {
    if (options_.ic3_reset_interval_
        && num_check_sat_since_reset_ >= options_.ic3_reset_interval_) {
      reset_solver();
    }

    const ProofGoal * pg = get_top_proof_goal();
    if(pg->idx > fidx)
      return true;

    if (is_blocked(pg)) {
      logger.log(3,
                 "Skipping already blocked proof goal <{}, {}>",
                 pg->target.term->to_string(),
                 pg->idx);
      remove_top_proof_goal();
      continue;
    };

    // block can fail, which just means a
    // new proof goal will be added in block
    if (block(pg)) {
      // if successfully blocked, then remove that proof goal
      // expecting the top proof goal to still be pg
      assert(pg == get_top_proof_goal());
      remove_top_proof_goal();
    } else if (!pg->idx) {
      // if a proof goal cannot be blocked at zero
      // then there's a counterexample
      // NOTE: creating a new allocation
      //       because the pg memory is already managed
      //       by proof_goals_
      
      // pop all added goals
      pg = get_top_proof_goal();
      while(true) {
        remove_top_proof_goal();
        if (has_proof_goals()) {
          pg = get_top_proof_goal();
          if (pg->idx > fidx)
            break;
          else
            continue;
        } else
          break;        
      }
      return false;
    }
  }
  assert(!has_proof_goals());
  return true;
} // block_all_goal_at_or_before

// ---------------------------------------------------------------------

void SygusPdr::disable_all_labels() {
  #define NOT(x) (solver_->make_term(Not, (x)))
  solver_->assert_formula(NOT(init_label_));
  solver_->assert_formula(NOT(trans_label_));
  for(const auto & fl : frame_labels_)
    solver_->assert_formula(NOT(fl));
}

// called in step to produce a bad state
// this is actually F /\ T -> bad'
bool SygusPdr::intersects_bad()
{
  push_solver_context();
  // assert the last frame (conjunction over clauses)
  assert_frame_labels(reached_k_ + 1);
  assert_trans_label();
  // see if it intersects with bad
  solver_->assert_formula(bad_next_);
  Result r = check_sat();

  if (r.is_sat()) {
    const IC3Formula &c = generalize_predecessor(
      reached_k_ + 1, IC3Formula(bad_, {bad_}, false) );
    // reduce c
    add_proof_goal(c, reached_k_ + 1, NULL);
  }

  pop_solver_context();

  assert(!r.is_unknown());
  return r.is_sat();
}

ProverResult SygusPdr::step_0()
{
  logger.log(1, "Checking if initial states satisfy property");
  assert(reached_k_ < 0);

  push_solver_context(); // sat(Init /\ bad)?
  solver_->assert_formula(init_label_);
  solver_->assert_formula(bad_);
  Result r = check_sat();
  if (r.is_sat()) {
    IC3Formula c = get_initial_bad_model();
    cex_pg_ = new ProofGoal(c, 0, nullptr);
    pop_solver_context();
    return ProverResult::FALSE;
  } else {
    assert(r.is_unsat());
    reached_k_ = 0;  // keep reached_k_ aligned with number of frames
  }
  pop_solver_context();

  push_solver_context(); // sat(Init /\ T /\ bad')?
  solver_->assert_formula(init_label_);
  assert_trans_label();
  solver_->assert_formula( bad_next_ );
  Result r = check_sat();
  if (r.is_sat()) {
    IC3Formula c = generalize_predecessor(0, IC3Formula(bad_, {bad_}, false));
    cex_pg_ = new ProofGoal(c, 0, nullptr);
    pop_solver_context();
    return ProverResult::FALSE;
  } else {
    assert(r.is_unsat());
    reached_k_ = 0;  // keep reached_k_ aligned with number of frames
  }
  pop_solver_context();

  return ProverResult::UNKNOWN;
} // step_0


// this is the same as the base class
// but because step_0 is not virtual
// I have to write it again
ProverResult SygusPdr::step(int i)
{
  if (i <= reached_k_) {
    return ProverResult::UNKNOWN;
  }

  if (reached_k_ < 0) {
    return SygusPdr::step_0();
  }

  // reached_k_ is the number of transitions that have been checked
  // at this point there are reached_k_ + 1 frames that don't
  // intersect bad, and reached_k_ + 2 frames overall
  assert(reached_k_ + 2 == frames_.size());
  logger.log(1, "Blocking phase at frame {}", i);
  // blocking phase
  while (intersects_bad()) {
    assert(has_proof_goals());
    if (!block_all()) {
      // counter-example
      return ProverResult::FALSE;
    }
  }

  logger.log(1, "Propagation phase at frame {}", i);
  // propagation phase
  push_frame();
  // in propagation, we should already have this "not(bad)" propagated
  for (size_t j = 1; j < frames_.size() - 1; ++j) {
    if (propagate(j)) {
      assert(j + 1 < frames_.size());
      // save the invariant
      // which is the frame that just had all terms
      // from the previous frames propagated
      invar_ = get_frame_term(j + 1);
      return ProverResult::TRUE;
    }
  }

  ++reached_k_;

  return ProverResult::UNKNOWN;
} // step


// this is the same as the base class
// but because step is not virtual
// I have to write it again
ProverResult SygusPdr::check_until(int k)
{
  initialize();
  // make sure derived class implemented initialize and called
  // this version of initialize with super::initialize or
  // (for experts only) set the initialized_ flag without
  // ever initializing base classes
  assert(initialized_);

  ProverResult res;
  RefineResult ref_res;
  int i = reached_k_ + 1;
  assert(i >= 0);
  while (i <= k) {
    // reset cex_pg_ to null
    // there might be multiple abstract traces if there's a derived class
    // doing abstraction refinement
    if (cex_pg_) {
      delete cex_pg_;
      cex_pg_ = nullptr;
    }

    res = SygusPdr::step(i);
    ref_res = REFINE_NONE;  // just a default value

    if (res == ProverResult::TRUE) {
      return res;
    } else if (res == ProverResult::FALSE) {
      // expecting cex_pg_ to be non-null and point to the first proof goal in a
      // trace
      assert(cex_pg_->target.term);
      ref_res = refine();
      if (ref_res == RefineResult::REFINE_NONE) {
        // found a concrete counterexample
        return res;
      } else if (ref_res == RefineResult::REFINE_FAIL) {
        logger.log(1, "Failed in refinement.");
        return ProverResult::UNKNOWN;
      }
    }

    // two cases
    // got unknown, so keep going
    // got false, but was able to refine successfully
    assert(res == ProverResult::UNKNOWN
           || (res == ProverResult::FALSE
               && ref_res == RefineResult::REFINE_SUCCESS));

    // increment i, unless there was a refinement step just done
    if (ref_res != RefineResult::REFINE_SUCCESS) {
      i++;
    }
  }

  return ProverResult::UNKNOWN;
}

}  // namespace pono


