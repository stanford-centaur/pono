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

#include <cassert>

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


void reduce_unsat_core_to_fixedpoint(
  const smt::Term & formula, 
  const smt::TermVec & assumptions,
  smt::TermVec & out, const smt::SmtSolver & reducer_) {
  // already pushed outside (because we want to disable all labels)
  reducer_->assert_formula(formula);

  // exit if the formula is unsat without assumptions.
  smt::Result r = reducer_->check_sat();
  if (r.is_unsat())
    return;

  out = (assumptions);
  while(true) {
    r = reducer_->check_sat_assuming(out);
    assert(r.is_unsat());

    smt::UnorderedTermSet core_set;
    reducer_->get_unsat_core(core_set);
    if (core_set.size() == out.size()) {
      break; // fixed point is reached
    }
    assert(core_set.size() < out.size());
    out.clear();
    out.insert(out.end(), core_set.begin(), core_set.end());    
  }
}

void reduce_unsat_core_linear(
    const smt::Term & formula,
    const smt::TermVec & assumptions,
    smt::TermVec & out,
    const smt::SmtSolver & reducer_) {

  // already pushed outside (because we want to disable all labels)
  reducer_->assert_formula(formula);

  // exit if the formula is unsat without assumptions.
  smt::Result r = reducer_->check_sat();
  if (r.is_unsat())
    return;

  r = reducer_->check_sat_assuming(assumptions);
  assert(r.is_unsat());
  size_t assump_pos_for_removal = 0;

  out = assumptions;
  while(assump_pos_for_removal < out.size()) {
    smt::TermVec assump_for_query;
    assump_for_query.reserve(out.size()-1);
    for (size_t idx = 0; idx < out.size(); ++ idx) {
      if (idx != assump_pos_for_removal)
        assump_for_query.push_back(out.at(idx));
    }

    r = reducer_->check_sat_assuming(assump_for_query);
    if (r.is_sat()) {
      ++ assump_pos_for_removal;
    } else {
      smt::UnorderedTermSet core_set;
      reducer_->get_unsat_core(core_set);
      { // remove those not in core_set from out
        // and we want to keep the previous order
        smt::TermVec new_assump;
        new_assump.reserve(core_set.size());
        for (const auto & l : out) {
          if (core_set.find(l) != core_set.end())
            new_assump.push_back(l);
        }
        out.swap(new_assump); // do "out = new_assump" w.o. copy
      }
      assert(!out.empty());
      assert(out.size() <= assump_for_query.size());
      // we don't need to change assump_pos_for_removal, because the size of
      // out is reduced by at least one, so in the next round
      // it is another element sitting at this location
    }
  } // end of while
} // end of reduce_unsat_core_linear


}  // namespace syntax_analysis
}  // namespace pono
