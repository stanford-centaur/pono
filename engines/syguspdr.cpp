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

using namespace smt;
using namespace std;

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

unsigned SygusPdr::GetScore(const smt::Term & t) {
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
  for (const auto & partial2full : to_full_model_map_) {
    if (partial2full.second)
      delete partial2full.second; // delete the full models
  }
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
    new CustomFunctionalTransitionSystem(orig_ts_, smt::TermTranslator( orig_ts_.solver()) );
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
  to_next_func_ = [&] (const smt::Term & v) -> smt::Term {
    return custom_ts_->next(v);
  };

  { // now create TermLearner
    term_learner_.reset(new syntax_analysis::TermLearner(
      custom_ts_->trans(), to_next_func_, solver_, 
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
    smt::TermVec * out_inputs,
    smt::TermVec * out_nexts) const {
  assert(false); // should not use this
  throw PonoException(
            "SyGuS PDR must be used with generalize_predecessor.");
  return IC3Formula();
} // SygusPdr::get_model_ic3formula


std::vector<IC3Formula> SygusPdr::inductive_generalization(
    size_t i, const IC3Formula & c) {
  // this is to produce the lemma
  // find the model
  auto model_pos = model2cube_.find( c.term );
  assert(model_pos != model2cube_.end());
  syntax_analysis::PerCexInfo & per_cex_info = setup_cex_info(*model_pos);
  // find the related info
  
  syntax_analysis::PredConstructor pred_collector_(
    to_next_func_, reducer_, per_cex_info
  );

  bool insufficient_pred = false;
  do {
    { // step 1 - check if the pred are good if not construct the model
      const auto & pred_nxt_ = pred_collector_.GetAllPredNext();
      push_solver_context();
      // F[i-1]
      assert_frame_labels(i - 1);
      // -c
      solver_->assert_formula(solver_->make_term(Not, c.term));
      // Trans
      assert_trans_label();
      // c'
      // solver_->assert_formula(ts_->next(c.term));
      for (const auto & pred_ : pred_nxt_) {
        solver_->assert_formula(ts_->next(pred_));
      }
      Result r = check_sat();
      if (r.is_sat()) {
        // TODO: extract model and do recursive block
        insufficient_pred = true;
        #error extractmodel
      }
      pop_solver_context();
    }

    if (insufficient_pred) {
      // recursive_block
    }
  }while(insufficient_pred);
  // at this point we have enough preds



  // will invoke block goal less than;

} // inductive_generalization

IC3Formula SygusPdr::generalize_predecessor(size_t i, const IC3Formula & c) {
  // used in rel_ind_check

  std::unordered_set<smt::Term> varlist;
  smt::Term bad_state_no_nxt = custom_ts_->to_next_func(c.term);
  if (has_assumptions) {
    constraints_curr_var_.push_back(bad_state_no_nxt);
    partial_model_getter.GetVarListForAsts(constraints_curr_var_, varlist);
    constraints_curr_var_.pop_back();
  } else {
    partial_model_getter.GetVarList(bad_state_no_nxt, varlist);
  }
  // extract the model based on var list
  
  // no need to pop
  // return the model and build IC3FormulaModel

  { // extract using keep_var_in_partial_model
    smt::Term conj_partial;
    smt::TermVec conjvec_partial;
    syntax_analysis::IC3FormulaModel::cube_t cube_partial;

    smt::Term conj_full;
    // smt::TermVec conjvec_full;
    syntax_analysis::IC3FormulaModel::cube_t cube_full;
    
    for (smt::Term v : varlist) {
      smt::Term val = solver_->get_value(v);
      auto eq = solver_->make_term(smt::Op(smt::PrimOp::Equal), v,val );
      { // full model part
        cube_full.emplace(v,val);
        // conjvec_full.push_back( eq );
        if (conj_full) {
          conj_full = solver_->make_term(smt::Op(smt::PrimOp::And), conj_full, eq);
        } else {
          conj_full = eq;
        }
      } // end of full model part
      if (keep_var_in_partial_model(v)) {
        cube_partial.emplace(v,val);
        conjvec_partial.push_back( eq );
        if (conj_partial) {
          conj_partial = solver_->make_term(smt::Op(smt::PrimOp::And), conj_partial, eq);
        } else {
          conj_partial = eq;
        }
      } // end of partial model
    }
    syntax_analysis::IC3FormulaModel * partial_model = 
      new syntax_analysis::IC3FormulaModel(std::move(cube_partial), conj_partial);
    syntax_analysis::IC3FormulaModel * full_model = 
      new syntax_analysis::IC3FormulaModel(std::move(cube_full), conj_full);

    assert(partial_model && full_model);
    model2cube_.emplace(conj_partial, partial_model);
    to_full_model_map_.emplace(partial_model, full_model);

    return IC3Formula(conj_partial, conjvec_partial, false /*not a disjunction*/ );\
  }
} // generalize_predecessor


bool SygusPdr::keep_var_in_partial_model(const smt::Term & v) const {
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




}  // namespace pono


