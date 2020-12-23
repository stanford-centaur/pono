/*********************                                                        */
/*! \file syntax_analysis.h
 ** \verbatim
 ** Top contributors (to current version):
 **   Hongce Zhang
 ** This file is part of the pono project.
 ** Copyright (c) 2020 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief Various functionalities for syntax analysis
 **
 ** 
 **/
 
#pragma once

#include "utils/syntax_analysis_common.h"


namespace pono {
namespace syntax_analysis {


class TermLearner;

class VarTermManager{
public:
  // type definition
  typedef std::unordered_map<std::string, PerVarsetInfo> terms_cache_t;
public:
  VarTermManager() {}
  void RegisterTermsToWalk(const smt::Term & t) { terms_to_check_.push_back(t); } // register the init by var and also trans by var, and property
  
  
  // this includes Constant Terms (will be inserted)
  const PerVarsetInfo & GetAllTermsForVarsInModel(Model * m, smt::SmtSolver & s);
  unsigned GetMoreTerms(Model * pre, Model * post, TermLearner & term_learner); // return delta terms

protected:
  std::unordered_set<std::string> constants_strings_;
  std::map<unsigned, std::vector<smt::Term>> width_to_constants_;  
  terms_cache_t terms_cache_; // include constants here
  
  std::vector<smt::Term> terms_to_check_;
  // from vars strings to the terms

protected:
  // 1. ----------------------------------------------
  // helps with the insertions
  void insert_from_constmap(const PerVarsetInfo::width_term_map_t & w_c_map) ;

  void term_const_w1_const(PerVarsetInfo & term_cache_item /*OUT*/, smt::SmtSolver & solver_);
  void const_to_per_varset(PerVarsetInfo & term_cache_item /*OUT*/, 
    unsigned width_bound_low /*IN*/, unsigned width_bound_high /*IN*/, unsigned & nterm_walked);

  unsigned insert_from_termsmap_w_width(
    const std::map<unsigned, smt::TermVec> & terms /*IN*/, PerVarsetInfo & term_cache_item /*OUT*/ , 
    unsigned width_bound_low /*IN*/, unsigned width_bound_high /*IN*/) ;

  void insert_vars_and_extracts(PerVarsetInfo & term_cache_item /*OUT*/,
    const smt::UnorderedTermSet & varset, smt::SmtSolver & solver_  );

  void insert_vars_only(
    PerVarsetInfo & term_cache_item /*OUT*/ , const smt::UnorderedTermSet & varset);

  void insert_split(
    PerVarsetInfo & term_cache_item /*OUT*/ , const smt::UnorderedTermSet & varset,
    smt::SmtSolver & s) ;

  // 2. ----------------------------------------------
  // helps with the Terms
  const PerVarsetInfo & SetupTermsForVarModelNormal(
    Model * m, const std::string & canonical_string);

  const PerVarsetInfo & SetupTermsForVarModeExt(
    Model * m, const std::string & canonical_string, smt::SmtSolver & s);

  const PerVarsetInfo & SetupTermsForVarModelVC(
    Model * m, const std::string & canonical_string);
  
  const PerVarsetInfo & SetupTermsForVarModeSplit(
    Model * m, const std::string & canonical_string, smt::SmtSolver & solver_);
  
}; // class VarTermManager

} // namespace syntax_analysis
} // namespace pono

