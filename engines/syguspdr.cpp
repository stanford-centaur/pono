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

    next_inputvars_.insert(next_inputvar);
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

SygusPdr::SygusPdr(Property & p, SolverEnum se) : super(p, se),
  partial_model_getter(solver_), custom_ts_(NULL)
{
  solver_->set_opt("produce-unsat-cores", "true");
}

SygusPdr::SygusPdr(Property & p, const SmtSolver & slv)
    : super(p, slv),
      partial_model_getter(solver_),
      custom_ts_(NULL)
{
  solver_->set_opt("produce-unsat-cores", "true");
}

SygusPdr::SygusPdr(const PonoOptions & opt,
                             Property & p,
                             const SolverEnum se)
    : super(opt, p, se),
      partial_model_getter(solver_),
      custom_ts_(NULL)
{
  solver_->set_opt("produce-unsat-cores", "true");
}

SygusPdr::SygusPdr(const PonoOptions & opt,
                             Property & p,
                             const SmtSolver & slv)
    : super(opt, p, slv),
      partial_model_getter(solver_),
      custom_ts_(NULL)
{
  solver_->set_opt("produce-unsat-cores", "true");
}

// ----------------------------------------------------------------
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
  
  // clear the mass if used by previous object
  // syntax_analysis::TermLearner::ClearCache();
  // partial to full model
  

  { // 1. register terms to find exprs
    // 2. extract parent from the same terms
    for (auto && v_nxtexpr_pair : custom_ts_->state_updates()) {
      sygus_term_manager_.RegisterTermsToWalk(v_nxtexpr_pair.second);
      parent_of_terms_.WalkBFS(v_nxtexpr_pair.second);
    }
    sygus_term_manager_.RegisterTermsToWalk(custom_ts_->init());
    parent_of_terms_.WalkBFS(custom_ts_->init());

    for (const auto & c : custom_ts_->constraints()) {
      sygus_term_manager_.RegisterTermsToWalk(c);
      parent_of_terms_.WalkBFS(c);
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



}  // namespace pono


