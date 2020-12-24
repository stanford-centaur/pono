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


// init 
SygusPdr::~SygusPdr() {
  for (const auto & ic3f_mptr : model2cube_) {
    if(ic3f_mptr.second)
      delete ic3f_mptr.second;
  }
} // SygusPdr::~SygusPdr

void SygusPdr::initialize()
{
  
  if (initialized_)
    return;

  super::initialize();
  
  // initialize the caches  
  // extract the operators
  op_extract_ = std::make_unique<OpExtractor>();
  op_extract_->WalkBFS(ts_msat_.init());
  op_extract_->WalkBFS(ts_msat_.trans());
  op_extract_->GetSyntaxConstruct().RemoveConcat();
  op_extract_->GetSyntaxConstruct().RemoveExtract();
  op_extract_->GetSyntaxConstruct().AndOrConvert();
  op_extract_->GetSyntaxConstruct().RemoveUnusedStructure();
  
  // clear the mass if used by previous object
  unsat_enum::Enumerator::ClearCache();
  unsat_enum::ParentExtract::ClearCache();
  unsat_enum::TermLearner::ClearCache();
  unsat_enum::TermScore::ClearCache();
  
  // clear the mass if used by previous object
  unsat_enum::Enumerator::ClearCache();
  unsat_enum::ParentExtract::ClearCache();
  unsat_enum::TermLearner::ClearCache();
  unsat_enum::TermScore::ClearCache();

  { // 1. register terms to find exprs
    // 2. extract parent from the same terms
    unsat_enum::ParentExtract parent_relation_extractor;
    for (auto && v_nxtexpr_pair : ts_.state_updates()) {
      sygus_term_manager_.RegisterTermsToWalk(v_nxtexpr_pair.second);
      parent_relation_extractor.WalkBFS(v_nxtexpr_pair.second);
    }
    sygus_term_manager_.RegisterTermsToWalk(ts_.init());
    parent_relation_extractor.WalkBFS(ts_.init());

    sygus_term_manager_.RegisterTermsToWalk(ts_.constraint());
    parent_relation_extractor.WalkBFS(ts_.constraint());

    sygus_term_manager_.RegisterTermsToWalk(property_.prop());
    parent_relation_extractor.WalkBFS(property_.prop());
  }

  { // now create TermLearner
    term_learner_.reset(new unsat_enum::TermLearner(
      ts_.trans(), to_next_func_, solver_, 
      unsat_enum::ParentExtract::GetParentRelation()
    ));
  }
}


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
}

}  // namespace pono


