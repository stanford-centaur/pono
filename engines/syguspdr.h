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
#include "modifiers/op_abstractor.h"
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

  // curr_var -> replace by var
  smt::Term to_next_func(const smt::Term & term) { 
      return FunctionalTransitionSystem::to_next_func(term); }
  const smt::UnorderedTermMap & input_var_to_next_map() const {return next_inputvars_;}

 protected:
  // next variables for the system inputs 
  smt::UnorderedTermMap next_inputvars_;

}; // class CustomFunctionalTransitionSystem

class SygusPdr : public IC3Base
{
 public:

  SygusPdr(const Property & p, const TransitionSystem & ts,
                const smt::SmtSolver & s,
                PonoOptions opt = PonoOptions());
  virtual ~SygusPdr(); // need to free the pointers

  // we cannot allow copy assignment
  SygusPdr & operator=(const SygusPdr &) = delete;

  typedef IC3Base super;
 

  // -----------------------------------------------------------------
  // override existing facilities to achieve bad = model(F /\ T -> P')
  // -----------------------------------------------------------------

 public:
  // some override that were unnecessary but I found that I need to do
  // inorder to adjust the style
  virtual ProverResult check_until(int k) override;

 protected:
  /** Perform a IC3 step
   *  @param i
   */
  ProverResult step(int i); // will be called in the parent version of check_until

  /** Perform the base IC3 step (zero case)
   */
  ProverResult step_0();  // will be called in the parent version of check_until

  // -----------------------------------------------------------------
  // Below are for May-Block
  // -----------------------------------------------------------------
  
  bool rel_ind_check_may_block(size_t i,
                              const IC3Formula & c,
                              IC3Formula & out);

  bool try_recursive_block_goal(const IC3Formula & to_block, unsigned fidx);

  // -----------------------------------------------------------------
  // Not to make a new model if blockable
  // -----------------------------------------------------------------

  bool rel_ind_check(size_t i,
                     const IC3Formula & c,
                     IC3Formula & out,
                     bool get_pred);
  
  bool block_all();

  // -----------------------------------------------------------------
  // pure virtual method implementations
  // -----------------------------------------------------------------

  virtual IC3Formula get_model_ic3formula() const override;

  bool ic3formula_check_valid(const IC3Formula & u) const override;


  virtual IC3Formula inductive_generalization(size_t i, const IC3Formula & c) override;

  virtual void predecessor_generalization(size_t i,
                                          const IC3Formula & c,
                                          IC3Formula & pred) override;

  void check_ts() const override;

  bool intersects_bad(IC3Formula & out) override;

  void initialize() override;


  virtual void abstract() override;

  virtual RefineResult refine() override;

  // -----------------------------------------------------------------
  // SyGuS related functions
  // -----------------------------------------------------------------

  // store the relation between IC3Formula to IC3Models
  std::unordered_map<smt::Term, syntax_analysis::IC3FormulaModel *> model2cube_;
  //std::unordered_map<syntax_analysis::IC3FormulaModel *, syntax_analysis::IC3FormulaModel *>
  //  to_full_model_map_;
  syntax_analysis::cex_term_map_t cex_term_map_;

  PartialModelGen partial_model_getter;

  IC3Formula get_initial_bad_model();
  std::pair<IC3Formula, syntax_analysis::IC3FormulaModel *>
    ExtractPartialModel(const smt::Term & p);
  std::pair<IC3Formula, syntax_analysis::IC3FormulaModel *>
    ExtractInitPrimeModel(const smt::Term & p_prime);
  std::tuple<IC3Formula, syntax_analysis::IC3FormulaModel *, syntax_analysis::IC3FormulaModel *>
    ExtractPartialAndFullModel(syntax_analysis::IC3FormulaModel * post_model);

  std::unique_ptr<syntax_analysis::OpExtractor> op_extract_;

  syntax_analysis::TermScore term_score_walker_;
  unsigned GetScore(const smt::Term & t);

  syntax_analysis::ParentExtract parent_of_terms_;
  syntax_analysis::VarTermManager sygus_term_manager_;
  std::unique_ptr<syntax_analysis::TermLearner> term_learner_;
  syntax_analysis::to_next_t to_next_func_;
  CustomFunctionalTransitionSystem * custom_ts_;
  
  smt::Term bad_next_;
  smt::TermVec constraints_curr_var_;

  bool has_assumptions;
  bool keep_var_in_partial_model(const smt::Term & v) const;
  void disable_all_labels();
  bool propose_new_terms(
    syntax_analysis::IC3FormulaModel * pre_model,
    syntax_analysis::IC3FormulaModel * post_model,
    const smt::Term & F_T_not_cex,
    const smt::Term & Init_prime,
    bool failed_at_init);

  syntax_analysis::PerCexInfo & setup_cex_info (syntax_analysis::IC3FormulaModel * post_model);
  IC3Formula select_predicates(const smt::Term & base, const smt::TermVec & preds_nxt);
  
  std::unique_ptr<OpAbstractor> op_abstractor_;
}; // class SygusPdr

}  // namespace pono

