/*********************                                                  */
/*! \file partial_model.cpp
** \verbatim
** Top contributors (to current version):
**   Hongce Zhang
** This file is part of the pono project.
** Copyright (c) 2020 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Class for performing dynamic cone of influence reduction based
**        on the model from the solver. This is essentially extracting a
**        partial model.
**
**/

#include "utils/partial_model.h"
#include "utils/str_util.h"

#include "assert.h"

namespace pono {

IC3Formula PartialModelGen::GetPartialModel(const smt::Term & ast) {
  GetVarList(ast);

  smt::Term conj;
  smt::TermVec conjvec;
  for (smt::Term v : dfs_vars_) {
    smt::Term val = solver_->get_value(v);
    auto eq = solver_->make_term(smt::Op(smt::PrimOp::Equal), v,val );
    conjvec.push_back( eq );
    if (conj) {
      conj = solver_->make_term(smt::Op(smt::PrimOp::And), conj, eq);
    } else {
      conj = eq;
    }
  }
  return IC3Formula(conj, conjvec, false /*not a disjunction*/ );
} // GetPartialModel to ic3formula


std::pair<IC3Formula,syntax_analysis::IC3FormulaModel> 
    PartialModelGen::GetPartialModelInCube(const smt::Term & ast) {
  
  GetVarList(ast);

  std::unordered_map<smt::Term, smt::Term> cube;

  smt::Term conj;
  smt::TermVec conjvec;
  for (smt::Term v : dfs_vars_) {
    smt::Term val = solver_->get_value(v);
    cube.emplace(v,val);
    auto eq = solver_->make_term(smt::Op(smt::PrimOp::Equal), v,val );
    conjvec.push_back( eq );
    if (conj) {
      conj = solver_->make_term(smt::Op(smt::PrimOp::And), conj, eq);
    } else {
      conj = eq;
    }
  }

  return std::make_pair(IC3Formula(conj, conjvec,
      false /*not a disjunction*/ ), 
    syntax_analysis::IC3FormulaModel(std::move(cube), conj));
}

void PartialModelGen::GetVarList(const smt::Term & ast ) {
  dfs_walked_.clear();
  dfs_vars_.clear();
  dfs_walk(ast);
}


void PartialModelGen::GetVarList(const smt::Term & ast, 
  std::unordered_set<smt::Term> & out_vars ) {

  dfs_walked_.clear();
  dfs_vars_.clear();
  dfs_walk(ast);
  out_vars.insert(dfs_vars_.begin(), dfs_vars_.end());
}

void PartialModelGen::GetVarListForAsts(const smt::TermVec & asts, 
  smt::UnorderedTermSet & out_vars ) {
  dfs_walked_.clear();
  dfs_vars_.clear();
  for (const auto & ast : asts)
    dfs_walk(ast);
  out_vars.insert(dfs_vars_.begin(), dfs_vars_.end());
}


/* Internal Function */
bool static extract_decimal_width(const std::string & s,
  std::string & decimal, std::string & width)
{
  if(s.substr(0,5) != "(_ bv")
    return false;
  auto space_idx = s.find(' ', 5);
  if(space_idx == std::string::npos )
    return false;
  auto rpara_idx = s.find(')', space_idx);
  if(rpara_idx == std::string::npos )
    return false;

  decimal = s.substr(5,space_idx-5);
  width = s.substr(space_idx+1,rpara_idx-(space_idx+1));
  assert(!width.empty());
  return true;
}


static std::string get_all_one(unsigned width) {
  std::vector<char> out = {1};

  for (unsigned idx = 1; idx < width; ++idx) {
    syntax_analysis::mul2(out);
    syntax_analysis::add1(out);
  }

  std::string ret;
  for (auto pos = out.rbegin(); pos != out.rend(); ++pos) {
    ret += *pos + '0';
  }
  return ret;
}

bool static convert_to_boolean_and_check(
  const std::string & decimal, const std::string & width, bool _0or1) {

  static std::unordered_map<unsigned, std::string> width2fullones;

  if (!_0or1) {
    for (auto c : decimal)
      if (c != '0')
        return false;
    return true;
  } // if 0, requires 0,00,000

  auto width_i = syntax_analysis::StrToULongLong(width,10);
  assert(width_i != 0);
  auto fullone_pos = width2fullones.find(width_i);
  if (fullone_pos == width2fullones.end()) {
    fullone_pos = width2fullones.emplace( width_i , get_all_one(width_i) ).first;
  }
  return decimal == fullone_pos->second;
}

/* Internal Function */
static inline bool is_all_zero(const std::string & s)  {
  if (s == "true")
    return false;
  if (s == "false")
    return true;
  assert (s.length() > 2);
  if (s.substr(0,2) == "#b") {
    for (auto pos = s.begin()+2; pos != s.end(); ++ pos)
      if (*pos != '0')
        return false;
    return true;
  } // else
  std::string decimal, width;
  bool conv_succ = extract_decimal_width(s, decimal, width);
  assert(conv_succ);

  return convert_to_boolean_and_check(decimal, width, false);
}

/* Internal Function */
static inline bool is_all_one(const std::string & s, uint64_t w)  {
  if (s == "true")
    return true;
  if (s == "false")
    return false;

  assert (s.length() > 2);
  if (s.substr(0,2) == "#b") {
    assert (s.length() - 2 <= w);
    if (s.length() - 2 < w) // if it has fewer zeros
      return false;
    for (auto pos = s.begin()+2; pos != s.end(); ++ pos)
      if (*pos != '1')
        return false;
    return true;  
  }
  std::string decimal, width;
  bool conv_succ = extract_decimal_width(s, decimal, width);
  assert(conv_succ);

  return convert_to_boolean_and_check(decimal, width, true);
}

/* Internal Macros */
#define ARG2(a1,a2)            \
      auto ptr = ast->begin();    \
      auto a1  = *(ptr++);      \
      auto a2  = *(ptr++);      \
      assert (ptr == ast->end()); 

#define ARG3(a1,a2,a3)            \
      auto ptr = ast->begin();    \
      auto a1  = *(ptr++);      \
      auto a2  = *(ptr++);      \
      auto a3  = *(ptr++);      \
      assert (ptr == ast->end());


void PartialModelGen::dfs_walk(const smt::Term & input_ast ) {
  smt::TermVec node_stack_;
  node_stack_.push_back(input_ast);
  while(!node_stack_.empty()) {
    const auto & ast = node_stack_.back();
    if (dfs_walked_.find(ast) != dfs_walked_.end()) {
      node_stack_.pop_back();
      continue;
    }
    dfs_walked_.insert(ast);

    smt::Op op = ast->get_op();
    if (op.is_null()) { // this is the root node
      if (ast->is_symbolic_const()) {
        dfs_vars_.insert(ast);
      }
      node_stack_.pop_back(); // no need to wait for the next time
      continue;
    } else { // non variable/non constant case
      if (op.prim_op == smt::PrimOp::Ite)  {
        ARG3(cond, texpr, fexpr)
        auto cond_val = solver_->get_value(cond);
        assert(cond_val->is_value());
        if ( is_all_one(cond_val->to_string(),1) ) {
          node_stack_.push_back(cond);
          node_stack_.push_back(texpr);
        }
        else {
          node_stack_.push_back(cond);
          node_stack_.push_back(fexpr);
        }
      } else if (op.prim_op == smt::PrimOp::Implies) {
        ARG2(left,right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        if (!( is_all_one(cond_left->to_string(),1) )) // if it is false
          node_stack_.push_back(left);
        else if ( is_all_one(cond_right->to_string(), 1) ) {
          node_stack_.push_back(right);
        } else {
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }
      } else if (op.prim_op == smt::PrimOp::And) {
        ARG2(left,right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        if (!( is_all_one(cond_left->to_string(),1) )) // if it is false
          node_stack_.push_back(left);
        else if (!(is_all_one(cond_right->to_string(), 1))) {
          node_stack_.push_back(right);
        } else {
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }
      } else if (op.prim_op == smt::PrimOp::Or) {
        ARG2(left,right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        if (is_all_one(cond_left->to_string(),1)) // if it is true
          node_stack_.push_back(left);
        else if (is_all_one(cond_right->to_string(), 1)) {
          node_stack_.push_back(right);
        } else  { // it is 0, so both matter
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }
      } else if (op.prim_op == smt::PrimOp::BVAnd || op.prim_op == smt::PrimOp::BVNand) {
        ARG2(left,right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        std::string left_val = cond_left->to_string();
        std::string right_val = cond_right->to_string();

        if (is_all_zero(left_val)) // if all zeros
          node_stack_.push_back(left);
        else if (is_all_zero(right_val)) {
          node_stack_.push_back(right);
        } else { // it is 0, so both matter
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }

      } else if (op.prim_op == smt::PrimOp::BVOr  || op.prim_op == smt::PrimOp::BVNor) {
        ARG2(left,right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        std::string left_val = cond_left->to_string();
        std::string right_val = cond_right->to_string();

        if (is_all_one(left_val, left->get_sort()->get_width())) // if all ones
          node_stack_.push_back(left);
        else if (is_all_one(right_val, right->get_sort()->get_width())) {
          node_stack_.push_back(right);
        } else { // it is 0, so both matter
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }
      } else if (op.prim_op == smt::PrimOp::BVMul) {
        ARG2(left,right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        std::string left_val = cond_left->to_string();
        std::string right_val = cond_right->to_string();

        if (is_all_zero(left_val)) // if all zeros
          node_stack_.push_back(left);
        else if (is_all_zero(right_val)) {
          node_stack_.push_back(right);
        } else { // it is not 0, so both matter
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }
      }
      else {
        for (const auto & arg : *ast)
          node_stack_.push_back(arg);
      }
    } // end non-variable case
  } // while ( not empty )
} // end of PartialModelGen::dfs_walk


} // namespace pono
