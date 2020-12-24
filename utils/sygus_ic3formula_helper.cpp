/*********************                                                        */
/*! \file syguspdr_ic3formula_helper.cpp
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

#include "utils/sygus_ic3formula_helper.h"
#include "utils/str_util.h"


namespace pono {
namespace syntax_analysis {

// intentionally create linker error
// TODO: add this function in the future
// IC3FormulaModel::IC3FormulaModel(const IC3Formula & f) : expr_(f.term) {
//   // extract the thing
//   for (const auto & c : f.children) {
//     const auto & op = c->get_op();
//     Op.prim_op == smt::PrimOp(Equal)
//   }
// }

IC3FormulaModel::IC3FormulaModel(IC3FormulaModel && f) :
  cube_(std::move(f.cube_)), expr_(std::move(f.expr_)) {
}
  
IC3FormulaModel & IC3FormulaModel::operator=(IC3FormulaModel && other) {
  if (this != &other) {
    cube_ = std::move(other.cube_);
    expr_ = std::move(other.expr_);
  }
  return *this;
}

std::string IC3FormulaModel::vars_to_canonical_string() const {
  std::vector<std::string> vars;
  for (auto && v_val : cube_) {
    vars.push_back(v_val.first->to_string());
  }
  std::sort(vars.begin(), vars.end());
  return Join(vars, "?<*>?"); // hope it won't appear in the the var names
}

std::string IC3FormulaModel::to_string() const {
  if (cube_.empty())
    return "true";
  
  std::string ret;
  for (auto && v_val : cube_ ) {
    ret += v_val.first->to_string()  + "= " +
       v_val.second->to_string() ;
    ret += " , ";
  }
  return ret;
}

void IC3FormulaModel::get_varset(std::unordered_set<smt::Term> & varset) const {
  for (auto && v_val : cube_) {
    varset.insert(v_val.first);
  }
}



}  // namespace syntax_analysis
}  // namespace pono
