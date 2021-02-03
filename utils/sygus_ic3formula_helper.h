/*********************                                                        */
/*! \file syguspdr_ic3formula_helper.h
** \verbatim
** Top contributors (to current version):
**   Hongce Zhang
** This file is part of the pono project.
** Copyright (c) 2020 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Utility functions to get variables from ic3formula
**        and to print it
**
**/

#pragma once

#include "smt-switch/smt.h"

#include <string>

namespace pono {

namespace syntax_analysis {

// I really want to cache some of the results
// instead of deciding the variables again and again
// and only class of partial model can construct

class IC3FormulaModel {
 public:
  typedef std::unordered_map<smt::Term, smt::Term> cube_t;
  // because we need to assign the expr_ and cube_ 
  // I intentionally make this a rvalue reference
  // because very often we need to construct the map
  // first and want to avoid a copy
  IC3FormulaModel(cube_t && cube, const smt::Term & expr) : 
    cube_(cube), expr_(expr) {}
  
  // here we really do the extraction
  // IC3FormulaModel(const IC3Formula & f);
  
  // move constructor
  IC3FormulaModel(IC3FormulaModel && f);
    
  IC3FormulaModel & operator=(IC3FormulaModel && other);
  
 protected:
  cube_t cube_;
  smt::Term expr_;
  
 public:
  std::string vars_to_canonical_string() const;
  std::string vars_val_to_canonical_string() const;
  std::string to_string() const;
  void get_varset(std::unordered_set<smt::Term> & varset) const;
  smt::Term to_expr() const { return expr_; }
  
}; // IC3FormulaModel

// will be used in replacement of unsatcore reducer
void reduce_unsat_core_to_fixedpoint(const smt::Term & formula, smt::UnorderedTermSet & core_inout, const smt::SmtSolver & solver_);

void reduce_unsat_core_linear(const smt::Term & formula, smt::TermList & assumption_list, const smt::SmtSolver & solver_);


}  // namespace syntax_analysis
}  // namespace pono


