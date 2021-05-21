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

std::string IC3FormulaModel::vars_val_to_canonical_string() const {
  std::vector<std::pair<std::string,std::string>> vars_vals;
  for (auto && v_val : cube_) {
    vars_vals.push_back(std::make_pair(v_val.first->to_string(), v_val.second->to_string()));
  }
  std::sort(vars_vals.begin(), vars_vals.end());
  std::string ret;
  for (const auto & v_val_pair :  vars_vals)
    ret += v_val_pair.first + "?<*>?" + v_val_pair.second;
  return ret; // hope it won't appear in the the var names
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
  smt::UnorderedTermSet & core_inout,
  const smt::SmtSolver & reducer_) {
  // already pushed outside (because we want to disable all labels)
  reducer_->assert_formula(formula);

  // exit if the formula is unsat without assumptions.
  smt::Result r = reducer_->check_sat();
  if (r.is_unsat())
    return;

  while(true) {
    r = reducer_->check_sat_assuming_set(core_inout);
    assert(r.is_unsat());

    smt::UnorderedTermSet core_out;
    reducer_->get_unsat_assumptions(core_out);
    if (core_inout.size() == core_out.size()) {
      break; // fixed point is reached
    }
    assert(core_out.size() < core_inout.size());
    core_inout.swap(core_out);  // namely, core_inout = core_out,  but no need to copy
  }
} // reduce_unsat_core_to_fixedpoint

// a helper function

void remove_and_move_to_next(smt::TermList & pred_set, smt::TermList::iterator & pred_pos,
  const smt::UnorderedTermSet & unsatcore) {

  auto pred_iter = pred_set.begin(); // pred_pos;
  auto pred_pos_new = pred_set.begin();

  bool reached = false;
  bool next_pos_found = false;

  while( pred_iter != pred_set.end() ) {

    if (pred_iter == pred_pos) {
      assert (!reached);
      reached = true;
    }
    
    if (unsatcore.find(*pred_iter) == unsatcore.end()) {
      assert (reached);
      pred_iter = pred_set.erase(pred_iter);
    } else {
      if (reached && ! next_pos_found) {
        pred_pos_new = pred_iter;
        next_pos_found = true;
      }
      ++ pred_iter;
    }
  } // end of while

  assert(reached);
  if (! next_pos_found) {
    assert (pred_iter == pred_set.end());
    pred_pos_new = pred_iter;
  }
  pred_pos = pred_pos_new;
} // remove_and_move_to_next

void reduce_unsat_core_linear(
    const smt::Term & formula,
    smt::TermList & assumption_list,
    const smt::SmtSolver & reducer_) {
  
  // already pushed outside (because we want to disable all labels)
  reducer_->assert_formula(formula);

  // exit if the formula is unsat without assumptions.
  smt::Result r = reducer_->check_sat();
  if (r.is_unsat())
    return;

  r = reducer_->check_sat_assuming_list(assumption_list);
  assert(r.is_unsat());
  auto to_remove_pos = assumption_list.begin();

  while(to_remove_pos != assumption_list.end()) {
    smt::Term term_to_remove = *to_remove_pos;
    auto pos_after = assumption_list.erase(to_remove_pos);
    r = reducer_->check_sat_assuming_list(assumption_list);
    to_remove_pos = assumption_list.insert(pos_after, term_to_remove);
    
    if (r.is_sat()) {
      ++ to_remove_pos;
    } else { // if unsat, we can remove
      smt::UnorderedTermSet core_set;
      reducer_->get_unsat_assumptions(core_set);
      // below function will update assumption_list and to_remove_pos
      remove_and_move_to_next(assumption_list, to_remove_pos, core_set);
    }
  } // end of while
} // end of reduce_unsat_core_linear


}  // namespace syntax_analysis
}  // namespace pono
