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
#include "utils/container_shortcut.h"
#include "utils/logger.h"
#include "utils/sygus_predicate_constructor.h"

using namespace smt;

namespace pono {

//#define DEBUG
#ifdef DEBUG
  #define D(...) logger.log( __VA_ARGS__ )
  #define INFO(...) D(0, __VA_ARGS__)
#else
  #define D(...) do {} while (0)
  #define INFO(...) logger.log(3, __VA_ARGS__)
#endif

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
SygusPdr::SygusPdr(const Property & p, const TransitionSystem & ts,
        const SmtSolver & s,
        PonoOptions opt)
  : super(p, ts, s, opt),
    partial_model_getter(solver_),
    custom_ts_(NULL),
    has_assumptions(true) // most conservative way
{
  solver_->set_opt("produce-unsat-cores", "true");

  // we need to have the reset-assertion capability 
  if(solver_->get_solver_enum() == SolverEnum::BTOR)
    solver_->set_opt("base-context-1", "true");
  // actually the above is already done in 
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
  if (!options_.ic3_pregen_)
    throw PonoException("SyGuS PDR requires IC3 predecessor generalization");
  if (!options_.ic3_indgen_)
    throw PonoException("SyGuS PDR requires IC3 inductive generalization");
  
  if (initialized_)
    return;

  ts_ = orig_ts_;

  // I really need the prime variable for inputs
  // otherwise the corner cases are hard to handle...


  super::initialize(); // I don't need the trans->prop thing

  bad_next_ = ts_.next(bad_); // bad is only available after parent's init
  { // add P to F[1]
    auto prop = smart_not(bad_);
    SygusPdr::constrain_frame(1, IC3Formula(prop, {prop}, true), true);
    // this may be redundant given the base ic3 use trans_label-> prop
  }


  if (options_.sygus_use_operator_abstraction_) {
    if ( options_.sygus_use_operator_abstraction_ == 2 )
      op_abstractor_ = std::make_unique<OpUfAbstractor>(
        orig_ts_, ts_, std::unordered_set<PrimOp>({BVMul, BVUdiv, BVSdiv, BVSmod, BVSrem, BVUrem}));
    else
      op_abstractor_ = std::make_unique<OpInpAbstractor>(
        orig_ts_, ts_, std::unordered_set<PrimOp>({BVMul, BVUdiv, BVSdiv, BVSmod, BVSrem, BVUrem}),
        bad_, 0);
    // op_abstractor_ = std::make_unique<OpUfAbstractor>(
    //   orig_ts_, ts_, std::unordered_set<PrimOp>({BVMul, BVUdiv, BVSdiv, BVSmod, BVSrem, BVUrem}), bad_, 4);
    if (options_.sygus_term_mode_ == SyGuSTermMode::TERM_MODE_AUTO)
      options_.sygus_term_mode_ = op_abstractor_->has_abstracted() ? 
        (SyGuSTermMode::SPLIT_FROM_DESIGN) : (SyGuSTermMode::FROM_DESIGN_LEARN_EXT);
  }
  if (options_.sygus_term_mode_ == SyGuSTermMode::TERM_MODE_AUTO)
    options_.sygus_term_mode_ = SyGuSTermMode::FROM_DESIGN_LEARN_EXT;

  SygusPdr::reset_solver(); // reset trans if abstracted
  // and even if not, I don't need the trans->prop thing


  { // create custom_ts_
    TermTranslator tmp_translator( ts_.solver());
    custom_ts_ = 
      new CustomFunctionalTransitionSystem(ts_, tmp_translator );
    custom_ts_->make_nextvar_for_inputs();
  }

  // has_assumption
  has_assumptions = false;
  for (const auto & c_initnext : custom_ts_->constraints()) {
    if (!c_initnext.second)
      continue;
    has_assumptions = true;
    assert(custom_ts_->no_next(c_initnext.first));
    // if (no_next) {
    constraints_curr_var_.push_back(c_initnext.first);
    // translate input_var to next input_var
    // but the state var ...
    // we will get to next anyway
    constraints_curr_var_.push_back(
      custom_ts_->to_next_func(
        solver_->substitute(c_initnext.first,
          custom_ts_->input_var_to_next_map())));
    // } // else skip
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

    for (const auto & c_next_init : custom_ts_->constraints()) {
      if (!c_next_init.second)
        continue;
      assert(custom_ts_->no_next(c_next_init.first));
      //if (custom_ts_->no_next(c_next_init.first)) {
      sygus_term_manager_.RegisterTermsToWalk(c_next_init.first);
      parent_of_terms_.WalkBFS(c_next_init.first);
      //  }
    }

    // sygus_term_manager_.RegisterTermsToWalk(property_.prop());
    // parent_of_terms_.WalkBFS(property_.prop());
  }

  // cache two lambda functions for sygus enum
  to_next_func_ = [&] (const Term & v) -> Term {
    return custom_ts_->next(v);
  };

  score_func_ = [this] (const Term & v) -> unsigned {
    return this->GetScore(v);
  };

  { // now create TermLearner
    term_learner_.reset(new syntax_analysis::TermLearner(
      to_next_func_, score_func_,
      solver_,  // you need to make trans dynamic
      parent_of_terms_
    ));
  }

} // initialize


void SygusPdr::check_ts() const
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
  //check each state has a next function
  for (const auto & sv : ts_.statevars() ) {
    const auto & update = ts_.state_updates();
    if (update.find(sv) == update.end())
      throw PonoException("State var `" + sv->to_string() +"` has no next function assigned."
        " Will mess up SyGuS PDR's internal logic." );
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


IC3Formula SygusPdr::get_model_ic3formula() const {
  // originally used in intersects_bad, step_0, rel_ind_check
  // now only in rel_ind_check
  // but we require 
  return IC3Formula();
} // SygusPdr::get_model_ic3formula


IC3Formula SygusPdr::inductive_generalization(
    size_t i, const IC3Formula & c) {
  assert(i>0);
  assert(!c.disjunction); // we need it to be a model

  // this is to produce the lemma
  // find the model
  auto model_pos = model2cube_.find( c.term );
  assert(model_pos != model2cube_.end());
  syntax_analysis::IC3FormulaModel * post_model = model_pos->second;

  syntax_analysis::PredConstructor pred_collector_(
    to_next_func_, solver_,
    post_model, // IC3FormulaModel * cex 
    setup_cex_info(post_model)
  );

  Term Init_prime = custom_ts_->next(ts_.init());
  const Term & cex = c.term;
  //Term not_cex = smart_not(cex);
  Term cex_prime = custom_ts_->next(cex);
  Term Fprev = get_frame_term(i - 1);
  // since we add the property as a lemma, so we should be okay here
  Term T = ts_.trans();
  Term F_T_not_cex = make_and( {Fprev, T} ); // , not_cex
  Term base = 
    solver_->make_term(Or,
      F_T_not_cex,
      Init_prime);

  IC3Formula pre_formula;
  syntax_analysis::IC3FormulaModel * pre_model = NULL;
  syntax_analysis::IC3FormulaModel * pre_full_model = NULL;
  bool failed_at_init;

  bool insufficient_pred = false;
  do {
    // refresh the frame lemmas
    Fprev = get_frame_term(i - 1);
    T = ts_.trans();
    F_T_not_cex = make_and( {Fprev, T} ); // , not_cex
    base = 
      solver_->make_term(Or,
        F_T_not_cex,
        Init_prime);

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
        failed_at_init = // if model(init_prime) == 1
          !( syntax_analysis::eval_val( solver_->get_value(Init_prime)->to_string() ) < syntax_analysis::eval_val("#b1") );
        if (pre_full_model) {
          delete pre_full_model;
          pre_full_model = NULL;
        }
        if (!failed_at_init) {
          // the idea of full model is we want to keep the assignment to
          // some inputs that could be thrown away
          std::tie(pre_formula, pre_model, pre_full_model) =
            ExtractPartialAndFullModel( post_model );
          D(3, "Generated MAY-block model {}", pre_model->to_string());
        } else {
          // only pre_model and failed_at_init are used in later
          std::tie(pre_formula, pre_model) =
            ExtractInitPrimeModel( Init_prime ); // this can be only the prime vars
          D(3, "Generated MAY-block model (init) {}", pre_model->to_string());
        }
        // note Here you must give a formula with current variables
      } else
        insufficient_pred = false;
      pop_solver_context();
    } // end of step 1

    if (insufficient_pred) {
      // extract Proof goal from partial_model
      if (!failed_at_init) {
        D(3, "Try recursive block the above model.");
        if(try_recursive_block_goal(pre_formula, i-1))
          continue; // see if we have next that we may need to block
      } // if we failed at init, then we will anyway need to 
      // propose new terms
      D(3, "Cannot block MAY-block model {}", pre_model->to_string());
      D(3, "Back to F[{}]->[{}], get new terms.", i-1,i);

      if (!failed_at_init)
        assert(pre_full_model);
      else
        assert(!pre_full_model);

      bool succ = 
        propose_new_terms(failed_at_init ? pre_model : pre_full_model, 
          post_model,
          F_T_not_cex, Init_prime, failed_at_init);
      assert(succ);
      delete pre_full_model;
      pre_full_model = NULL;

      // it has loop inside
      // and it must succeed because eventually we will end with bit-level things
    }
  }while(insufficient_pred);

  if (pre_full_model)
    delete pre_full_model;

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
  D(3,"[IterativeReduction] Initial core size {}", preds_nxt.size());

  UnorderedTermSet pred_set(preds_nxt.begin(), preds_nxt.end());
  TermVec deduplicate_preds_nxt(pred_set.begin(), pred_set.end());

  if(solver_->get_solver_enum() == SolverEnum::BTOR) {
    push_solver_context();
    disable_all_labels();
    syntax_analysis::reduce_unsat_core_to_fixedpoint(base, deduplicate_preds_nxt, unsatcore, solver_);
    pop_solver_context();
  } else
    reducer_.reduce_assump_unsatcore(base, deduplicate_preds_nxt, unsatcore, NULL, 0, options_.random_seed_);

  D(3,"[IterativeReduction] End core size {}", unsatcore.size());
  std::vector<std::pair<unsigned, Term>> score_vec;
  for (const auto & t : unsatcore) {
    auto score = GetScore(t);
    score_vec.push_back(std::make_pair(score, t));
  }
  std::sort(score_vec.begin(),score_vec.end());
  std::unordered_map<Term, unsigned> term_id_map;
  unsigned pidx = 0;
  for (auto pos = score_vec.rbegin(); pos != score_vec.rend(); ++ pos) {
    sorted_unsatcore.push_back(pos->second);
    D(4, "{}: Score: {} expr: {}", pidx, pos->first, pos->second->to_string());
    term_id_map.emplace(pos->second, pidx ++);
  }
  unsatcore.clear(); // reuse this
  if(solver_->get_solver_enum() == SolverEnum::BTOR) {
    push_solver_context();
    disable_all_labels();
    syntax_analysis::reduce_unsat_core_linear(base, sorted_unsatcore, unsatcore, solver_);
    pop_solver_context();
    #ifdef DEBUG
      std::cout << "{";
      for (const auto & p : unsatcore) {
        auto pos = term_id_map.find(p); assert(pos != term_id_map.end());
        std::cout << pos->second << ",";
      }
      std::cout << "}" << std::endl;
    #endif
  } else 
    reducer_.linear_reduce_assump_unsatcore(base, sorted_unsatcore, unsatcore, NULL, 0);
  
  assert(!unsatcore.empty());

  TermVec curr_lits;
  curr_lits.reserve(unsatcore.size());
  for (const auto &l : unsatcore)  {
    curr_lits.push_back(custom_ts_->curr(l)); // make sure we also convert input back
  }
  return ic3formula_negate(ic3formula_conjunction(curr_lits));
} // select_predicates

bool SygusPdr::propose_new_terms(
    syntax_analysis::IC3FormulaModel * pre_model,
    syntax_analysis::IC3FormulaModel * post_model,
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
        ts_.trans(),
        failed_at_init, options_.sygus_term_mode_);
    D(3, "[propose-new-term] Round {}. Get {} new terms.", proposing_new_terms_round, n_new_terms);
    if (n_new_terms != 0) {
      syntax_analysis::PredConstructor pred_collector_(
        to_next_func_, solver_,
        post_model, // IC3FormulaModel * cex 
        setup_cex_info(post_model) 
        // you need to set this up every time you use it
        // this is needed to eval new terms
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
      auto res = check_sat();
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
            sygus_term_manager_.GetAllTermsForVarsInModel(post_model, solver_,
              options_.sygus_term_mode_,
              options_.sygus_term_extract_depth_,
              options_.sygus_initial_term_width_,
              options_.sygus_initial_term_inc_,
              options_.sygus_accumulated_term_bound_ ) ) );
  } // if no terms, set up them first

  syntax_analysis::PerCexInfo & per_cex_info = cex_term_map_pos->second;
  unsigned nterm = 0;

  // then evalute the terms on the cex
  push_solver_context();
  disable_all_labels();

  solver_->assert_formula( post_model->to_expr() );
  auto res = check_sat();
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
  TermVec conjvec;
  for (const auto & v : ts_.statevars()) {
    Term val = solver_->get_value(v);
    auto eq = solver_->make_term(Op(PrimOp::Equal), v,val );
    conjvec.push_back(eq);

    if (conj_full)
      conj_full = solver_->make_term(Op(PrimOp::And), conj_full, eq);
    else
      conj_full = eq;
  }

  for (const auto & v : ts_.inputvars()) {
    Term val = solver_->get_value(v);
    auto eq = solver_->make_term(Op(PrimOp::Equal), v,val );
    conjvec.push_back(eq);
    
    if (conj_full)
      conj_full = solver_->make_term(Op(PrimOp::And), conj_full, eq);
    else
      conj_full = eq;
  }
  assert(conj_full);
  assert(!conjvec.empty());

  return IC3Formula(conj_full, conjvec, false /*not a disjunction*/ );
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
    Term curr_v = custom_ts_->curr(v);
    auto eq = solver_->make_term(Op(PrimOp::Equal), curr_v ,val );
    cube_partial.emplace(curr_v,val);
    conjvec_partial.push_back( eq );
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
    new syntax_analysis::IC3FormulaModel(std::move(cube_partial), conj_partial);

  assert(partial_model);
  model2cube_.emplace(conj_partial, partial_model);  

  return std::make_pair(
    IC3Formula(conj_partial, conjvec_partial, false /*not a disjunction*/ ),
    partial_model
  );
} // ExtractInitPrimeModel 

std::tuple<IC3Formula, syntax_analysis::IC3FormulaModel *, syntax_analysis::IC3FormulaModel *>
    SygusPdr::ExtractPartialAndFullModel(syntax_analysis::IC3FormulaModel * post_model) {
  // we get varset from post_model, this will save us from one term walk
  UnorderedTermSet varset;
  post_model->get_varset(varset);
  TermVec pre_var;
  for (const auto & v : varset) {
    if (IN(v,ts_.statevars())) {
      auto update_pos = ts_.state_updates().find(v);
      assert( update_pos != ts_.state_updates().end() );
      pre_var.push_back(update_pos->second);
    } // if it is input then we should not need
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
      auto eq = solver_->make_term(Op(PrimOp::Equal), v,val );

      cube_full.emplace(v, val);
      if (conj_full)
        conj_full = solver_->make_term(Op(PrimOp::And), conj_full, eq);
      else
        conj_full = eq;

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
    assert(conj_full); // at least there should be var in the model
    if (conj_full == nullptr)
      conj_full = solver_true_;

    if (conj_partial == nullptr) {
      conj_partial = solver_true_;
      assert(conjvec_partial.empty());
      conjvec_partial.push_back(solver_true_);
    }
    syntax_analysis::IC3FormulaModel * partial_model = 
      new syntax_analysis::IC3FormulaModel(std::move(cube_partial), conj_partial);

    syntax_analysis::IC3FormulaModel * full_model = 
      new syntax_analysis::IC3FormulaModel(std::move(cube_full), conj_full);

    assert(partial_model && full_model);
    model2cube_.emplace(conj_partial, partial_model);

    return std::make_tuple(
      IC3Formula(conj_partial, conjvec_partial, false /*not a disjunction*/ ),
      partial_model,
      full_model
    );
  } // end model making block
} // ExtractPartialAndFullModel

std::pair<IC3Formula, syntax_analysis::IC3FormulaModel *>
    SygusPdr::ExtractPartialModel(const Term & p) {
  // extract using keep_var_in_partial_model  
  assert(custom_ts_->no_next(p));

  UnorderedTermSet varlist;
  Term bad_state_no_nxt = custom_ts_->to_next_func(
    solver_->substitute(p, custom_ts_->input_var_to_next_map()));
  // we need to make sure input vars are mapped to next input vars

  D(4, "[PartialModel] prime state : {}", bad_state_no_nxt->to_string());
  if (has_assumptions) {
    D(4, "[PartialModel] assumptions (mapped): {}",constraints_curr_var_.size());
    unsigned idx = 0;
    for (const auto & c : constraints_curr_var_)
      D(4, "[PartialModel] assumption #{} : {}", idx ++, c->to_string());
    constraints_curr_var_.push_back(bad_state_no_nxt);
    partial_model_getter.GetVarListForAsts(constraints_curr_var_, varlist);
    constraints_curr_var_.pop_back();
  } else {
    partial_model_getter.GetVarList(bad_state_no_nxt, varlist);
  }

  {
    D(4, "[PartialModel] before cutting vars: ");
    for (const auto & v : varlist)
      D(4, "[PartialModel] {} := {} ", v->to_string(), solver_->get_value(v)->to_string());
    D(4, "[PartialModel] ------------------- ");
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
  if (conj_partial == nullptr) {
    conj_partial = solver_true_;
    assert(conjvec_partial.empty());
    conjvec_partial.push_back(solver_true_);
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


void SygusPdr::predecessor_generalization(size_t i, const Term & cterm, IC3Formula & pred) {
  // used in rel_ind_check
  // extract the model based on var list
  // NOTE: i may be incorrect, it is given as F/\T->(here is i)

  // no need to pop (pop in rel_ind_check)
  // return the model and build IC3FormulaModel
  auto partial_full_model = ExtractPartialModel(cterm);
  pred = partial_full_model.first;
} // generalize_predecessor


bool SygusPdr::keep_var_in_partial_model(const Term & v) const {
  if (options_.sygus_use_operator_abstraction_) {
    assert (op_abstractor_ != nullptr);
    const auto & dummy_inputs = op_abstractor_->dummy_inputs();
    if (dummy_inputs.find(v) != dummy_inputs.end())
      return true;
  }

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
  if (!options_.sygus_use_operator_abstraction_)
    return REFINE_NONE;

  assert(witness_length()>0);

  TermVec new_constraints;
  bool succ = op_abstractor_->refine_with_constraints(
    cex_, bad_, new_constraints);
  assert(succ == (!new_constraints.empty()));
  if (succ) {
    for (const auto & c : new_constraints) {
      ts_.add_constraint(c, true);
      custom_ts_->add_constraint(c, true);
      D(4, "[Refine] : {} ", c->to_string());

      // rework the constraints buffers
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
    reset_solver(); // because tr/init could be different now

    // rework the constraints
    has_assumptions = ! (custom_ts_->constraints().empty());

    return REFINE_SUCCESS;
  } // else
  return REFINE_NONE;
} // SygusPdr::refine()


// ------------------------------------------------------------------
//  May Block 

/**
 * Priority queue of proof obligations inspired by open-source ic3ia
 * implementation
 */
class MayProofGoalQueue
{
 public:
  ~MayProofGoalQueue() { clear(); }

  void clear()
  {
    for (auto p : store_) {
      delete p;
    }
    store_.clear();
    while (!queue_.empty()) {
      queue_.pop();
    }
  }

  void new_proof_goal(const IC3Formula & c,
                      unsigned int t,
                      const ProofGoal * n = NULL)
  {
    ProofGoal * pg = new ProofGoal(c, t, n);
    queue_.push(pg);
    store_.push_back(pg);
  }

  ProofGoal * top() { return queue_.top(); }
  void pop() { queue_.pop(); }
  bool empty() const { return queue_.empty(); }

 private:
  typedef std::
      priority_queue<ProofGoal *, std::vector<ProofGoal *>, ProofGoalOrder>
          Queue;
  Queue queue_;
  std::vector<ProofGoal *> store_;
};

bool SygusPdr::try_recursive_block_goal(const IC3Formula & to_block, unsigned fidx) {
  assert(!solver_context_);
  MayProofGoalQueue proof_goals;
  proof_goals.new_proof_goal(to_block, fidx, nullptr);
  while(!proof_goals.empty()) {
    const ProofGoal * pg = proof_goals.top();
    if (pg->idx == 0) { // fail
      // do we refine ? no, because this is just a may proof goal
      return false;
    }

    if (is_blocked(pg)) {
      assert(pg == proof_goals.top());
      proof_goals.pop();
      continue;
    }

    //  try and see if it is blockable
    IC3Formula collateral;  // populated by rel_ind_check: if unsat collateral:=target  if sat collateral:=partial model
    if (rel_ind_check_may_block(pg->idx, pg->target, collateral)) {
      // this proof goal can be blocked
      assert(!solver_context_);
      assert(collateral.term);
      logger.log(
          3, "Blocking term at frame {}: {}", pg->idx, pg->target.term);

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
      logger.log(
          3, "Blocking term at frame {}: {} , \n --- with {}", pg->idx, pg->target.term, collateral.term->to_string());
    } else {
      if (collateral.term == solver_true_)
        return false; // when we intersect with 0
      // could not block this proof goal
      assert(collateral.term);
      proof_goals.new_proof_goal(collateral, pg->idx - 1, pg);
    }
  } // end of while
  assert(proof_goals.empty());
  return true;
} // try_recursive_block_goal_at_or_before



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

  solver_->assert_formula(custom_ts_->next(c.term));

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
    out.children = {solver_true_};
    out.disjunction = false;
    pop_solver_context();
    return r.is_unsat();
  }
  
  out = c;
  assert(!out.disjunction); // we need it to be a model
  pop_solver_context();
  assert(!r.is_unknown());
  return r.is_unsat();
} // rel_ind_check_may_block

// ---------------------------------------------------------------------
//   Override for: bad = model(F/\T -> P')
// ---------------------------------------------------------------------

void SygusPdr::disable_all_labels() {
  #define NOT(x) (solver_->make_term(Not, (x)))
  solver_->assert_formula(NOT(init_label_));
  solver_->assert_formula(NOT(trans_label_));
  solver_->assert_formula(NOT(bad_label_));
  for(const auto & fl : frame_labels_)
    solver_->assert_formula(NOT(fl));
  #undef NOT
}

// called in step to produce a bad state
// this is actually F /\ T -> bad'
bool SygusPdr::reaches_bad(IC3Formula & out)
{
  push_solver_context();
  // assert the last frame (conjunction over clauses)
  assert_frame_labels(frontier_idx());
  assert_trans_label();
  // see if it intersects with bad
  solver_->assert_formula(bad_next_);
  Result r = check_sat();

  if (r.is_sat()) {
    predecessor_generalization(
      frames_.size(), bad_, out );
    // actually the first is not used
  }

  pop_solver_context();

  assert(!r.is_unknown());
  return r.is_sat();
} // intersects_bad


ProverResult SygusPdr::step_01()
{ // differs from base class now in get_initial_bad_model
  assert(reached_k_ < 1);
  if (reached_k_ < 0) {
    logger.log(1, "Checking if initial states satisfy property");

    push_solver_context();
    solver_->assert_formula(init_label_);
    solver_->assert_formula(bad_);
    Result r = check_sat();
    if (r.is_sat()) {
      pop_solver_context();
      // trace is only one bad state that intersects with initial
      cex_.clear();
      cex_.push_back(bad_);
      return ProverResult::FALSE;
    } else {
      assert(r.is_unsat());
      reached_k_ = 0;  // keep reached_k_ aligned with number of frames
    }
    pop_solver_context();
  }

  assert(reached_k_ == 0);
  logger.log(1, "Checking if property can be violated in one-step");

  push_solver_context();
  solver_->assert_formula(init_label_);
  solver_->assert_formula(trans_label_);
  solver_->assert_formula(custom_ts_->next(bad_));
  Result r = check_sat();
  if (r.is_sat()) {
    IC3Formula c = get_initial_bad_model();
    pop_solver_context();
    ProofGoal * pg = new ProofGoal(c, 0, nullptr);
    SygusPdr::reconstruct_trace(pg, cex_);
    delete pg;
    return ProverResult::FALSE;
  } else {
    assert(r.is_unsat());
    reached_k_ = 1;  // keep reached_k_ aligned with number of frames
  }
  pop_solver_context();

  return ProverResult::UNKNOWN;
}

// difference : 1. no need for invar to be conjuncted
// 2. call SygusPdr::block_all
ProverResult SygusPdr::step(int i)
{
  if (i <= reached_k_) {
    return ProverResult::UNKNOWN;
  }

  if (reached_k_ < 1) {
    return step_01();
  }

  // reached_k_ is the number of transitions that have been checked
  // at this point there are reached_k_ + 1 frames that don't
  // intersect bad, and reached_k_ + 2 frames overall
  assert(reached_k_ == frontier_idx());
  logger.log(1, "Blocking phase at frame {}", i);
  if (!SygusPdr::block_all()) {
    // counter-example
    return ProverResult::FALSE;
  }

  logger.log(1, "Propagation phase at frame {}", i);
  // propagation phase
  push_frame();
  for (size_t j = 1; j < frontier_idx(); ++j) {
    if (SygusPdr::propagate(j)) {
      assert(j + 1 < frames_.size());
      // save the invariant
      // which is the frame that just had all terms
      // from the previous frames propagated
      invar_ = get_frame_term(j + 1);
      // invar_ = solver_->make_term(And, invar_, smart_not(bad_));
      return ProverResult::TRUE;
    }
  }

  reset_solver();

  ++reached_k_;

  return ProverResult::UNKNOWN;
}


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
  assert(reached_k_ + 1 >= 0);
  while (i <= k) {
    res = step(i);

    if (res == ProverResult::FALSE) {
      assert(cex_.size());
      RefineResult s = refine();
      if (s == REFINE_SUCCESS) {
        continue;
      } else if (s == REFINE_NONE) {
        // this is a real counterexample
        assert(cex_.size());
        return ProverResult::FALSE;
      } else {
        assert(s == REFINE_FAIL);
        throw PonoException("Refinement failed");
      }
    } else {
      ++i;
    }

    if (res != ProverResult::UNKNOWN) {
      return res;
    }
  }

  return ProverResult::UNKNOWN;
}// check_until


// ---------------------------------------------------------------------
//   Override to not make a new model when unsat
//   otherwise the lemma generation will not get the ic3model
// ---------------------------------------------------------------------

// The difference is (1) not to create a new model when unsat
// (2) use custom_ts->next (which also has next var for input)
// (3) also check intersection with init
bool SygusPdr::rel_ind_check(size_t i,
                            const IC3Formula & c,
                            IC3Formula & out,
                            bool get_pred)
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
  solver_->assert_formula(solver_->make_term(Not, c.term));
  // Trans
  assert_trans_label();

  // use assumptions for c' so we can get cheap initial
  // generalization if the check is unsat
  solver_->assert_formula(custom_ts_->next(c.term));
  Result r = check_sat();
  if (r.is_sat()) {
    if (get_pred) {
      out = get_model_ic3formula();
      if (options_.ic3_pregen_) {
        predecessor_generalization(i, c.term, out);
        assert(out.term);
        assert(out.children.size());
        assert(!out.disjunction);  // expecting a conjunction
      }
    }
    assert(ic3formula_check_valid(out));
    pop_solver_context();
    if (get_pred) {
      assert(out.term);
      assert(out.children.size());

      // in case of abstraction, maybe not
      if (!options_.sygus_use_operator_abstraction_) {
        // this check needs to be here after the solver context has been popped
        // if i == 1 and there's a predecessor, then it should be an initial state
        assert(i != 1 || check_intersects_initial(out.term));

        // should never intersect with a frame before F[i-1]
        // otherwise, this predecessor should have been found
        // in a previous step (before a new frame was pushed)
        assert(i < 2 || !check_intersects(out.term, get_frame_term(i - 2)));
      }

    }
    return r.is_unsat();
  } // else
  pop_solver_context();
  push_solver_context();
  
    solver_->assert_formula(init_label_);
    solver_->assert_formula(c.term);
    r = check_sat();
    if (r.is_sat()) {
      // get model at 0
      out.term = solver_true_;
      out.children = {solver_true_};
      out.disjunction = false;
      pop_solver_context();
      return r.is_unsat();
    }

  pop_solver_context();
  // unsat
  out = c;
  
  assert(!solver_context_);
  assert(!r.is_unknown());
  return r.is_unsat();
} // rel_ind_check


// there is no difference from the base class
// but I want to update rel_ind_check
bool SygusPdr::block_all()
{
  assert(!solver_context_);
  MayProofGoalQueue proof_goals;
  IC3Formula goal;
  while (SygusPdr::reaches_bad(goal)) {
    assert(goal.term);            // expecting non-null
    assert(proof_goals.empty());  // bad should be the first goal each iteration
    proof_goals.new_proof_goal(goal, frontier_idx(), nullptr);

    while (!proof_goals.empty()) {
      const ProofGoal * pg = proof_goals.top();

      if (!pg->idx) {
        // went all the way back to initial
        // need to create a new proof goal that's not managed by the queue
        SygusPdr::reconstruct_trace(pg, cex_);

        // in case this is spurious, clear the queue of proof goals
        // which might not have been precise
        // TODO might have to change this if there's an algorithm
        // that refines but can keep proof goals around
        proof_goals.clear();

        return false;
      }

      if (is_blocked(pg)) {
        logger.log(3,
                   "Skipping already blocked proof goal <{}, {}>",
                   pg->target.term,
                   pg->idx);
        // remove the proof goal since it has already been blocked
        assert(pg == proof_goals.top());
        proof_goals.pop();
        continue;
      }

      IC3Formula collateral;  // populated by rel_ind_check
      if (rel_ind_check(pg->idx, pg->target, collateral, true)) {
        // this proof goal can be blocked
        assert(!solver_context_);
        assert(collateral.term);
        logger.log(
            3, "Blocking term at frame {}: {}", pg->idx, pg->target.term);

        // remove the proof goal now that it has been blocked
        assert(pg == proof_goals.top());
        proof_goals.pop();

        if (options_.ic3_indgen_) {
          collateral = inductive_generalization(pg->idx, collateral);
        } else {
          // just negate the term
          collateral = ic3formula_negate(collateral);
        }

        size_t idx = find_highest_frame(pg->idx, collateral);
        assert(idx >= pg->idx);

        assert(collateral.disjunction);
        assert(collateral.term);
        assert(collateral.children.size());
        SygusPdr::constrain_frame(idx, collateral);

        // re-add the proof goal at a higher frame if not blocked
        // up to the frontier
        if (idx < frontier_idx()) {
          assert(!pg->target.disjunction);
          proof_goals.new_proof_goal(pg->target, idx + 1, pg->next);
        }

      } else {
        // could not block this proof goal
        assert(collateral.term);
        if (collateral.term == solver_true_)
          ; // do nothing // because we need to extend it to 0 to 
        proof_goals.new_proof_goal(collateral, pg->idx - 1, pg);
      }
    }  // end while(!proof_goals.empty())

    assert(!(goal = IC3Formula()).term);  // in debug mode, reset it
  }                                       // end while(reaches_bad(goal))

  assert(proof_goals.empty());
  return true;
} // block_all

// -----------------------------------------------------------------
// allow input in the constraint (when there are additional constraints)
// -----------------------------------------------------------------


static bool subsumes(const IC3Formula &a, const IC3Formula &b)
{
  assert(a.disjunction);
  assert(a.disjunction == b.disjunction);
  const TermVec &ac = a.children;
  const TermVec &bc = b.children;
  // NOTE: IC3Formula children are sorted on construction
  //       Uses unique id of term, from term->get_id()
  return ac.size() <= bc.size()
         && std::includes(bc.begin(), bc.end(), ac.begin(), ac.end());
}


void SygusPdr::constrain_frame(size_t i, const IC3Formula & constraint,
                              bool new_constraint)
{
  assert(solver_context_ == 0);
  assert(i < frame_labels_.size());
  assert(constraint.disjunction);
  if (has_assumptions || options_.sygus_use_operator_abstraction_ == 1) {
    // options_.sygus_use_operator_abstraction_ == 2 will not introduce new inputs
    // that must be kept
    assert(ts_.no_next(constraint.term));
  } else
    assert(ts_.only_curr(constraint.term));

  if (new_constraint) {
    for (size_t j = 1; j <= i; ++j) {
      std::vector<IC3Formula> & Fj = frames_.at(j);
      size_t k = 0;
      for (size_t l = 0; l < Fj.size(); ++l) {
        if (!subsumes(constraint, Fj[l])) {
          Fj[k++] = Fj[l];
        }
      }
      Fj.resize(k);
    }
  }

  assert(i > 0);  // there's a special case for frame 0

  constrain_frame_label(i, constraint);
  frames_.at(i).push_back(constraint);
}

// -----------------------------------------------------------------
// Require the capability to reset solver
// -----------------------------------------------------------------

void SygusPdr::reset_solver()
{
  assert(solver_context_ == 0);

  // for syguspdr, we must be able to reset solver
  assert (!failed_to_reset_solver_);

  try {
    solver_->reset_assertions();

    // Now need to add back in constraints at context level 0
    logger.log(2, "SygusPdr: Reset solver and now re-adding constraints.");

    // define init, trans, and bad labels
    assert(init_label_ == frame_labels_.at(0));
    solver_->assert_formula(
        solver_->make_term(Implies, init_label_, ts_.init()));

    solver_->assert_formula(
        solver_->make_term(Implies, trans_label_, ts_.trans()));

    // assume property in pre-state

    solver_->assert_formula(solver_->make_term(Implies, bad_label_, bad_));

    for (size_t i = 0; i < frames_.size(); ++i) {
      for (const auto & constraint : frames_.at(i)) {
        constrain_frame_label(i, constraint);
      }
    }
  }
  catch (SmtException & e) {
    logger.log(1,
               "Failed to reset solver (underlying solver must not support "
               "it). Disabling solver resets for rest of run.");
    throw e; // SyGuS requires reset solver capability 
  }

  num_check_sat_since_reset_ = 0;
}

// difference : I don't need to conjuct with property, 
// because it is already added as a lemma/constraint
Term SygusPdr::get_frame_term(size_t i) const
{
  // TODO: decide if frames should hold IC3Formulas or terms
  //       need to special case initial state if using IC3Formulas
  if (i == 0) {
    // F[0] is always the initial states constraint
    return ts_.init();
  }

  Term res = solver_true_;
  for (size_t j = i; j < frames_.size(); ++j) {
    for (const auto &u : frames_[j]) {
      res = solver_->make_term(And, res, u.term);
    }
  }
  return res;
}

// Difference : I need it to call the new constrain_frame
// which allow input in the constraints
bool SygusPdr::propagate(size_t i)
{
  assert(!solver_context_);
  assert(i < frontier_idx());

  std::vector<IC3Formula> & Fi = frames_.at(i);

  size_t k = 0;
  IC3Formula gen;
  for (size_t j = 0; j < Fi.size(); ++j) {
    const IC3Formula & c = Fi.at(j);
    assert(c.disjunction);
    assert(c.term);
    assert(c.children.size());

    // NOTE: rel_ind_check works on conjunctions
    //       need to negate
    if (rel_ind_check(i + 1, ic3formula_negate(c), gen, false)) {
      // can push to next frame
      // got unsat-core based generalization
      assert(gen.term);
      assert(gen.children.size());
      SygusPdr::constrain_frame(i + 1, ic3formula_negate(gen), false);
    } else {
      // have to keep this one at this frame
      Fi[k++] = c;
    }
  }

  // get rid of garbage at end of frame
  Fi.resize(k);

  return Fi.empty();
} // propagate

// difference: lift the restriction on only_curr vars
void SygusPdr::reconstruct_trace(const ProofGoal * pg, TermVec & out)
{
  assert(!solver_context_);
  assert(pg);
  assert(pg->target.term);
  assert(check_intersects_initial(pg->target.term));

  out.clear();
  while (pg) {
    out.push_back(pg->target.term);
    if (!has_assumptions && !options_.sygus_use_operator_abstraction_)
      assert(ts_.only_curr(out.back()));
    else
      assert(ts_.no_next(out.back()));
    pg = pg->next;
  }

  // always add bad as last state so it's a full trace
  // NOTE this is because the reaches_bad implementation
  out.push_back(bad_);
}

}  // namespace pono


