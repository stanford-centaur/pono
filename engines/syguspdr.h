/*********************                                                  */
/*! \file syguspdr.h
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
**        based on
**
**        Syntax-Guided Synthesis for Lemma Generation in Unbounded Hardware Verification
**            -- Hongce Zhang, Aarti Gupta, Sharad Malik
**
**        and the open source implementation:
**
**        https://github.com/zhanghongce/cosa2/tree/sygus-apdr
**
**
**/

#pragma once

#include "engines/ic3base.h"
#include "core/fts.h"
#include "utils/sygus_ic3formula_helper.h"
#include "utils/partial_model.h"
#include "utils/syntax_analysis_common.h"
#include "utils/syntax_analysis_walker.h"
#include "utils/syntax_analysis.h"

namespace pono {

// this class is added simply because SygusPdr would like to see such a transition system
class CustomFunctionalTransitionSystem : public FunctionalTransitionSystem {
 public:
  CustomFunctionalTransitionSystem() : FunctionalTransitionSystem() { }
  CustomFunctionalTransitionSystem(const smt::SmtSolver & s) : FunctionalTransitionSystem(s) { }
  CustomFunctionalTransitionSystem(const TransitionSystem & other_ts,
                             smt::TermTranslator & tt) : FunctionalTransitionSystem(other_ts, tt) {  }
  void make_nextvar_for_inputs();
  smt::Term to_next_func(const smt::Term & term) { 
      return FunctionalTransitionSystem::to_next_func(term); }
 protected:
  // next variables for the system inputs 
  smt::UnorderedTermSet next_inputvars_;

}; // class CustomFunctionalTransitionSystem

class SygusPdr : public IC3Base
{
 public:
  SygusPdr(Property & p, smt::SolverEnum se);
  SygusPdr(Property & p, const smt::SmtSolver & s);
  SygusPdr(const PonoOptions & opt, Property & p, smt::SolverEnum se);
  SygusPdr(const PonoOptions & opt,
                Property & p,
                const smt::SmtSolver & s);
  virtual ~SygusPdr(); // need to free the pointers

  // we cannot allow copy assignment
  SygusPdr & operator=(const SygusPdr &) = delete;

  typedef IC3Base super;
 
 protected:

  // pure virtual method implementations

  IC3Formula get_model_ic3formula(
      smt::TermVec * out_inputs = nullptr,
      smt::TermVec * out_nexts = nullptr) const override;

  bool ic3formula_check_valid(const IC3Formula & u) const override;

  std::vector<IC3Formula> inductive_generalization(
      size_t i, const IC3Formula & c) override;

  IC3Formula generalize_predecessor(size_t i, const IC3Formula & c) override;

  void check_ts() const override;

  bool intersects_bad() override;

  void initialize() override;

  // store the relation between IC3Formula to IC3Models
  std::unordered_map<smt::Term, syntax_analysis::IC3FormulaModel *> model2cube_;
  std::unordered_map<syntax_analysis::IC3FormulaModel *, syntax_analysis::IC3FormulaModel *>
    to_full_model_map_;
  syntax_analysis::cex_term_map_t cex_term_map_;

  PartialModelGen partial_model_getter;
  std::unique_ptr<syntax_analysis::OpExtractor> op_extract_;

  syntax_analysis::TermScore term_score_walker_;
  unsigned GetScore(const smt::Term & t);

  syntax_analysis::ParentExtract parent_of_terms_;
  syntax_analysis::VarTermManager sygus_term_manager_;
  std::unique_ptr<syntax_analysis::TermLearner> term_learner_;
  syntax_analysis::to_next_t to_next_func_;
  
}; // class SygusPdr

}  // namespace pono

