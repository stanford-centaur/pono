/*********************                                                        */
/*! \file sygus_predicate_constructor.h
** \verbatim
** Top contributors (to current version):
**   Hongce Zhang
** This file is part of the pono project.
** Copyright (c) 2020 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Utility class to help manage/construct the predicates
**
**/

#pragma once

#include "utils/syntax_analysis.h"
#include "utils/sygus_ic3formula_helper.h"
#include "smt-switch/utils.h"

#include <vector>
#include <functional>
#include <list>


namespace pono {
namespace syntax_analysis {

class PredConstructor {
  
protected:
  smt::SmtSolver & solver_;
  // smt::Term trans_;
  // smt::Term init_;
  // smt::Term prev_;
  IC3FormulaModel * cex_; // the cexs to block
    
  
protected:
  // btor_var_to_msat_t btor_var_to_msat_;
  // const btor_var_to_msat_var_map_t & btor_var_to_msat_var_map_;
  to_next_t to_next_;
  
  PerCexInfo & per_cex_info_; // per var set info here
  
  void terms_to_predicates();
  
  void TermsDumping() const;

  smt::Term zero,one;

  smt::Term smart_EQ(const smt::Term &l, const  smt::Term & r);
  smt::Term smart_NEQ(const smt::Term &l, const  smt::Term & r);
  smt::Term smart_LT(const smt::Term &l, const  smt::Term & r);
  smt::Term smart_LE(const smt::Term &l, const  smt::Term & r);

public:
  // constraints are included in the T part
  PredConstructor(
    to_next_t to_next_func,
    smt::SmtSolver & solver,
    // const smt::Term & T_btor, const smt::Term & Init_btor, const smt::Term & Fprev_btor,
    IC3FormulaModel * cex,
    PerCexInfo & per_cex_info // from pdr class, this will allow us to not use static data
    //bool init_per_cex_info,
    // VarTermManager & var_term_extractor 
    // setup_cex_info will be in pdr class
  );
  
  const smt::TermVec & GetAllPredNext() const { return per_cex_info_.predicates_nxt; }

  uint64_t term_summary() const;

}; // PredConstructor

} // namespace syntax_analysis
} // namespace pono

