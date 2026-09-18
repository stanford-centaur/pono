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

#include <algorithm>
#include <cassert>
#include <list>

#include "utils/str_util.h"

namespace pono {

static inline uint64_t unified_width(const smt::Term & ast)
{
  return ast->get_sort()->get_sort_kind() == smt::SortKind::BOOL
             ? 1
             : ast->get_sort()->get_width();
}

IC3Formula PartialModelGen::GetPartialModel(const smt::Term & ast)
{
  GetVarList(ast);

  smt::Term conj;
  smt::TermVec conjvec;
  for (smt::Term v : dfs_vars_) {
    smt::Term val = solver_->get_value(v);
    auto eq = solver_->make_term(smt::Op(smt::PrimOp::Equal), v, val);
    conjvec.push_back(eq);
    if (conj) {
      conj = solver_->make_term(smt::Op(smt::PrimOp::And), conj, eq);
    } else {
      conj = eq;
    }
  }
  return IC3Formula(conj, conjvec, false /*not a disjunction*/);
}  // GetPartialModel to ic3formula

std::pair<IC3Formula, syntax_analysis::IC3FormulaModel>
PartialModelGen::GetPartialModelInCube(const smt::Term & ast)
{
  GetVarList(ast);

  std::unordered_map<smt::Term, smt::Term> cube;

  smt::Term conj;
  smt::TermVec conjvec;
  for (smt::Term v : dfs_vars_) {
    smt::Term val = solver_->get_value(v);
    cube.emplace(v, val);
    auto eq = solver_->make_term(smt::Op(smt::PrimOp::Equal), v, val);
    conjvec.push_back(eq);
    if (conj) {
      conj = solver_->make_term(smt::Op(smt::PrimOp::And), conj, eq);
    } else {
      conj = eq;
    }
  }

  return std::make_pair(
      IC3Formula(conj, conjvec, false /*not a disjunction*/),
      syntax_analysis::IC3FormulaModel(std::move(cube), conj));
}

IC3Formula PartialModelGen::GetPartialModel_bitlevel(const smt::Term & ast)
{
  auto high = static_cast<int>(unified_width(ast)) - 1;
  std::unordered_map<smt::Term, std::vector<std::pair<int, int>>>
      output_var_slices;

  GetVarListForAsts_in_bitlevel(
      { { ast, std::vector<std::pair<int, int>>(1, { high, 0 }) } },
      output_var_slices);

  smt::Term conj;
  smt::TermVec conjvec;
  for (const auto & var_vec_pair : output_var_slices) {
    const auto & v = var_vec_pair.first;
    const auto & slices = var_vec_pair.second;
    smt::Term val = solver_->get_value(v);
    if (v->get_sort()->get_sort_kind() == smt::SortKind::BOOL) {
      auto eq = solver_->make_term(smt::Op(smt::PrimOp::Equal), v, val);
      conjvec.push_back(eq);
      if (conj) {
        conj = solver_->make_term(smt::Op(smt::PrimOp::And), conj, eq);
      } else {
        conj = eq;
      }
      continue;
    }

    for (const auto & bitfield : slices) {
      auto left = solver_->make_term(
          smt::Op(smt::PrimOp::Extract, bitfield.first, bitfield.second), v);
      auto right = solver_->make_term(
          smt::Op(smt::PrimOp::Extract, bitfield.first, bitfield.second), val);
      auto eq = solver_->make_term(smt::Op(smt::PrimOp::Equal), left, right);

      conjvec.push_back(eq);
      if (conj) {
        conj = solver_->make_term(smt::Op(smt::PrimOp::And), conj, eq);
      } else {
        conj = eq;
      }
    }
  }
  if (!conj) {  // generally, this should not happen
    conj = solver_->make_term(true);
    conjvec.push_back(conj);
  }
  return IC3Formula(conj, conjvec, false /*not a disjunction*/);
}

std::pair<IC3Formula, syntax_analysis::IC3FormulaModel>
PartialModelGen::GetPartialModelInCube_bitlevel(const smt::Term & ast)
{
  auto high = static_cast<int>(unified_width(ast)) - 1;
  std::unordered_map<smt::Term, std::vector<std::pair<int, int>>>
      output_var_slices;

  GetVarListForAsts_in_bitlevel(
      { { ast, std::vector<std::pair<int, int>>(1, { high, 0 }) } },
      output_var_slices);

  std::unordered_map<smt::Term, smt::Term> cube;

  smt::Term conj;
  smt::TermVec conjvec;
  for (const auto & var_vec_pair : output_var_slices) {
    const auto & v = var_vec_pair.first;
    const auto & slices = var_vec_pair.second;
    smt::Term val = solver_->get_value(v);
    if (v->get_sort()->get_sort_kind() == smt::SortKind::BOOL) {
      auto eq = solver_->make_term(smt::Op(smt::PrimOp::Equal), v, val);
      cube.emplace(v, val);
      conjvec.push_back(eq);
      if (conj) {
        conj = solver_->make_term(smt::Op(smt::PrimOp::And), conj, eq);
      } else {
        conj = eq;
      }
      continue;
    }

    for (const auto & bitfield : slices) {
      auto left = solver_->make_term(
          smt::Op(smt::PrimOp::Extract, bitfield.first, bitfield.second), v);
      auto right = solver_->make_term(
          smt::Op(smt::PrimOp::Extract, bitfield.first, bitfield.second), val);
      auto eq = solver_->make_term(smt::Op(smt::PrimOp::Equal), left, right);

      cube.emplace(left, right);
      conjvec.push_back(eq);
      if (conj) {
        conj = solver_->make_term(smt::Op(smt::PrimOp::And), conj, eq);
      } else {
        conj = eq;
      }
    }
  }
  if (!conj) {
    conj = solver_->make_term(true);
    conjvec.push_back(conj);
  }

  return std::make_pair(
      IC3Formula(conj, conjvec, false /*not a disjunction*/),
      syntax_analysis::IC3FormulaModel(std::move(cube), conj));
}

void PartialModelGen::GetVarList(const smt::Term & ast)
{
  dfs_walked_.clear();
  dfs_vars_.clear();
  dfs_walk(ast);
}

void PartialModelGen::GetVarList(const smt::Term & ast,
                                 std::unordered_set<smt::Term> & out_vars)
{
  dfs_walked_.clear();
  dfs_vars_.clear();
  dfs_walk(ast);
  out_vars.insert(dfs_vars_.begin(), dfs_vars_.end());
}

static inline bool in_range(int i, int left, int right)
{
  return i <= left + 1 && i >= right - 1;
}

static inline bool mergeable(const std::pair<int, int> & a,
                             const std::pair<int, int> & b)
{
  return in_range(a.first, b.first, b.second)
         || in_range(a.second, b.first, b.second)
         || in_range(b.first, a.first, a.second)
         || in_range(b.second, a.first, a.second);
}

static inline std::pair<int, int> merge_range(
    const std::pair<int, int> & l, const std::pair<int, int> & r)
{
  return { std::max(l.first, r.first), std::min(l.second, r.second) };
}

static std::vector<std::pair<int, int>> merge_intervals(
    const std::vector<std::pair<int, int>> & intervals)
{
  if (intervals.size() <= 1) return intervals;

  std::list<std::pair<int, int>> merged;
  for (auto r : intervals) {
    auto pos = merged.begin();
    while (pos != merged.end()) {
      if (mergeable(*pos, r)) {
        r = merge_range(*pos, r);
        merged.erase(pos);
        pos = merged.begin();
      } else {
        ++pos;
      }
    }
    merged.push_back(r);
  }

  return std::vector<std::pair<int, int>>(merged.begin(), merged.end());
}

void PartialModelGen::GetVarListForAsts_in_bitlevel(
    const std::unordered_map<smt::Term, std::vector<std::pair<int, int>>> &
        input_asts_slices,
    std::unordered_map<smt::Term, std::vector<std::pair<int, int>>> & out)
{
  // Multiple paths can require overlapping or adjacent slices of the same
  // variable. Merge them so callers can build compact equality constraints.
  dfs_walked_extract_.clear();
  std::unordered_map<smt::Term, PairSet> varset_slices;
  for (const auto & ast_slice_pair : input_asts_slices) {
    for (const auto & h_l_pair : ast_slice_pair.second) {
      dfs_walk_bitlevel(
          ast_slice_pair.first, h_l_pair.first, h_l_pair.second, varset_slices);
    }
  }
  for (const auto & var_slice_pair : varset_slices) {
    std::vector<std::pair<int, int>> intervals(var_slice_pair.second.begin(),
                                               var_slice_pair.second.end());
    out.emplace(var_slice_pair.first, merge_intervals(intervals));
  }
}

void PartialModelGen::GetVarListForAsts(const smt::TermVec & asts,
                                        smt::UnorderedTermSet & out_vars)
{
  dfs_walked_.clear();
  dfs_vars_.clear();
  for (const auto & ast : asts) dfs_walk(ast);
  out_vars.insert(dfs_vars_.begin(), dfs_vars_.end());
}

/* Internal Function */
bool static extract_decimal_width(const std::string & s,
                                  std::string & decimal,
                                  std::string & width)
{
  if (s.substr(0, 5) != "(_ bv") return false;
  auto space_idx = s.find(' ', 5);
  if (space_idx == std::string::npos) return false;
  auto rpara_idx = s.find(')', space_idx);
  if (rpara_idx == std::string::npos) return false;

  decimal = s.substr(5, space_idx - 5);
  width = s.substr(space_idx + 1, rpara_idx - (space_idx + 1));
  assert(!width.empty());
  return true;
}

static std::string get_all_one(unsigned width)
{
  std::vector<char> out = { 1 };

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

bool static convert_to_boolean_and_check(const std::string & decimal,
                                         const std::string & width,
                                         bool _0or1)
{
  static std::unordered_map<unsigned, std::string> width2fullones;

  if (!_0or1) {
    for (auto c : decimal)
      if (c != '0') return false;
    return true;
  }  // if 0, requires 0,00,000

  auto width_i = syntax_analysis::StrToULongLong(width, 10);
  assert(width_i != 0);
  auto fullone_pos = width2fullones.find(width_i);
  if (fullone_pos == width2fullones.end()) {
    fullone_pos = width2fullones.emplace(width_i, get_all_one(width_i)).first;
  }
  return decimal == fullone_pos->second;
}

/* Internal Function */
static inline bool is_all_zero(const std::string & s)
{
  if (s == "true") return false;
  if (s == "false") return true;
  assert(s.length() > 2);
  if (s.substr(0, 2) == "#b") {
    for (auto pos = s.begin() + 2; pos != s.end(); ++pos)
      if (*pos != '0') return false;
    return true;
  }  // else
  std::string decimal, width;
  bool conv_succ = extract_decimal_width(s, decimal, width);
  assert(conv_succ);

  return convert_to_boolean_and_check(decimal, width, false);
}

/* Internal Function */
static inline bool is_all_one(const std::string & s, uint64_t w)
{
  if (s == "true") return true;
  if (s == "false") return false;

  assert(s.length() > 2);
  if (s.substr(0, 2) == "#b") {
    assert(s.length() - 2 <= w);
    if (s.length() - 2 < w)  // if it has fewer zeros
      return false;
    for (auto pos = s.begin() + 2; pos != s.end(); ++pos)
      if (*pos != '1') return false;
    return true;
  }
  std::string decimal, width;
  bool conv_succ = extract_decimal_width(s, decimal, width);
  assert(conv_succ);

  return convert_to_boolean_and_check(decimal, width, true);
}

static std::string decimal_to_binary(std::string decimal)
{
  // Avoid fixed-width integer conversion here; SMT models may print decimal
  // bit-vector constants wider than native integer types.
  decimal.erase(0, decimal.find_first_not_of('0'));
  if (decimal.empty()) return "0";

  std::string binary;
  while (!decimal.empty()) {
    std::string quotient;
    int carry = 0;
    for (char c : decimal) {
      int value = carry * 10 + (c - '0');
      int digit = value / 2;
      carry = value % 2;
      if (!quotient.empty() || digit != 0) {
        quotient.push_back(static_cast<char>('0' + digit));
      }
    }
    binary.push_back(static_cast<char>('0' + carry));
    decimal = quotient;
  }

  std::reverse(binary.begin(), binary.end());
  return binary;
}

// #b0000 -> 0000
// (_ bvX W) -> X in binary (W-wide)
static inline std::string to_unified_bvconst(const std::string & s)
{
  // Normalize solver-specific BV value strings so bit index 0 is the last
  // character before find_consecutive_zeros_ones reverses the string.
  if (s.substr(0, 2) == "#b") return s.substr(2);

  std::string decimal, width;
  bool conv_succ = extract_decimal_width(s, decimal, width);
  assert(conv_succ);

  auto binary = decimal_to_binary(decimal);
  auto width_int = syntax_analysis::StrToULongLong(width, 10);
  if (binary.length() < width_int) {
    binary.insert(0, width_int - binary.size(), '0');
  }
  return binary;
}

static void find_consecutive_zeros_ones(
    std::string s,
    const std::vector<std::pair<int, int>> & ranges,
    std::vector<std::pair<int, int>> & ranges_of_zeros,
    std::vector<std::pair<int, int>> & ranges_of_ones)
{
  std::reverse(s.begin(), s.end());
  for (const auto & msb_lsb_pair : ranges) {
    auto msb = msb_lsb_pair.first;
    auto lsb = msb_lsb_pair.second;
    assert(lsb >= 0 && msb < static_cast<int>(s.length()) && msb >= lsb);

    int prev_i = lsb;
    for (int i = lsb + 1; i <= msb + 1; ++i) {
      if (i == msb + 1 || s[i] != s[prev_i]) {
        if (s[prev_i] == '0') {
          ranges_of_zeros.push_back({ i - 1, prev_i });
        } else {
          ranges_of_ones.push_back({ i - 1, prev_i });
        }
        prev_i = i;
      }
    }
  }
}

static inline std::string bool2bv_str(const std::string & in)
{
  if (in == "true") return "#b1";
  if (in == "false") return "#b0";
  return in;
}

static int find_rightmost_index_with_different_bit(const std::string & left,
                                                   const std::string & right,
                                                   size_t width)
{
  auto l = to_unified_bvconst(bool2bv_str(left));
  auto r = to_unified_bvconst(bool2bv_str(right));
  assert(!l.empty() && !r.empty());
  assert(l.length() == r.length() && l.length() == width);

  int bitindex = 0;
  for (int pos = static_cast<int>(l.length()) - 1; pos >= 0; --pos) {
    if (l.at(pos) != r.at(pos)) return bitindex;
    ++bitindex;
  }
  assert(static_cast<size_t>(bitindex) <= width);
  return -1;
}

static int find_leftmost_index_with_different_bit(const std::string & left,
                                                  const std::string & right,
                                                  size_t width)
{
  auto l = to_unified_bvconst(bool2bv_str(left));
  auto r = to_unified_bvconst(bool2bv_str(right));
  assert(!l.empty() && !r.empty());
  assert(l.length() == r.length() && l.length() == width);

  int bitindex = static_cast<int>(width) - 1;
  for (size_t pos = 0; pos < width; ++pos) {
    if (l.at(pos) != r.at(pos)) return bitindex;
    --bitindex;
  }
  assert(bitindex >= -1);
  return 0;
}

/* Internal Macros */
#define ARG1(a1)             \
  auto ptr = ast->begin();   \
  auto a1 = *(ptr++);        \
  assert(ptr == ast->end());

#define ARG2(a1, a2)       \
  auto ptr = ast->begin(); \
  auto a1 = *(ptr++);      \
  auto a2 = *(ptr++);      \
  assert(ptr == ast->end());

#define ARG3(a1, a2, a3)   \
  auto ptr = ast->begin(); \
  auto a1 = *(ptr++);      \
  auto a2 = *(ptr++);      \
  auto a3 = *(ptr++);      \
  assert(ptr == ast->end());

void PartialModelGen::dfs_walk(const smt::Term & input_ast)
{
  smt::TermVec node_stack_;
  node_stack_.push_back(input_ast);
  while (!node_stack_.empty()) {
    const auto & ast = node_stack_.back();
    if (dfs_walked_.find(ast) != dfs_walked_.end()) {
      node_stack_.pop_back();
      continue;
    }
    dfs_walked_.insert(ast);

    smt::Op op = ast->get_op();
    if (op.is_null()) {  // this is the root node
      if (ast->is_symbolic_const()) {
        dfs_vars_.insert(ast);
      }
      node_stack_.pop_back();  // no need to wait for the next time
      continue;
    } else {  // non variable/non constant case
      if (op.prim_op == smt::PrimOp::Ite) {
        ARG3(cond, texpr, fexpr)
        auto cond_val = solver_->get_value(cond);
        assert(cond_val->is_value());
        if (is_all_one(cond_val->to_string(), 1)) {
          node_stack_.push_back(cond);
          node_stack_.push_back(texpr);
        } else {
          node_stack_.push_back(cond);
          node_stack_.push_back(fexpr);
        }
      } else if (op.prim_op == smt::PrimOp::Implies) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        if (!(is_all_one(cond_left->to_string(), 1)))  // if it is false
          node_stack_.push_back(left);
        else if (is_all_one(cond_right->to_string(), 1)) {
          node_stack_.push_back(right);
        } else {
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }
      } else if (op.prim_op == smt::PrimOp::And) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        if (!(is_all_one(cond_left->to_string(), 1)))  // if it is false
          node_stack_.push_back(left);
        else if (!(is_all_one(cond_right->to_string(), 1))) {
          node_stack_.push_back(right);
        } else {
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }
      } else if (op.prim_op == smt::PrimOp::Or) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        if (is_all_one(cond_left->to_string(), 1))  // if it is true
          node_stack_.push_back(left);
        else if (is_all_one(cond_right->to_string(), 1)) {
          node_stack_.push_back(right);
        } else {  // it is 0, so both matter
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }
      } else if (op.prim_op == smt::PrimOp::BVAnd
                 || op.prim_op == smt::PrimOp::BVNand) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        std::string left_val = cond_left->to_string();
        std::string right_val = cond_right->to_string();

        if (is_all_zero(left_val))  // if all zeros
          node_stack_.push_back(left);
        else if (is_all_zero(right_val)) {
          node_stack_.push_back(right);
        } else {  // it is 0, so both matter
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }

      } else if (op.prim_op == smt::PrimOp::BVOr
                 || op.prim_op == smt::PrimOp::BVNor) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        std::string left_val = cond_left->to_string();
        std::string right_val = cond_right->to_string();

        if (is_all_one(left_val, left->get_sort()->get_width()))  // if all ones
          node_stack_.push_back(left);
        else if (is_all_one(right_val, right->get_sort()->get_width())) {
          node_stack_.push_back(right);
        } else {  // it is 0, so both matter
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }
      } else if (op.prim_op == smt::PrimOp::BVMul) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        std::string left_val = cond_left->to_string();
        std::string right_val = cond_right->to_string();

        if (is_all_zero(left_val))  // if all zeros
          node_stack_.push_back(left);
        else if (is_all_zero(right_val)) {
          node_stack_.push_back(right);
        } else {  // it is not 0, so both matter
          node_stack_.push_back(left);
          node_stack_.push_back(right);
        }
      } else {
        for (const auto & arg : *ast) node_stack_.push_back(arg);
      }
    }  // end non-variable case
  }  // while ( not empty )
}  // end of PartialModelGen::dfs_walk

void PartialModelGen::dfs_walk_bitlevel(
    const smt::Term & input_ast,
    int high,
    int low,
    std::unordered_map<smt::Term, PairSet> & varset_slice)
{
  std::vector<std::pair<smt::Term, std::pair<int, int>>> node_stack_;

  // The stack stores a term plus the output slice whose value must be
  // justified by the partial model.
  auto input_width = static_cast<int>(unified_width(input_ast));
  assert(high >= low && low >= 0 && high < input_width);
  node_stack_.push_back({ input_ast, { high, low } });

  while (!node_stack_.empty()) {
    auto ast = node_stack_.back().first;
    auto extracted_bit = node_stack_.back().second;
    auto width = static_cast<int>(unified_width(ast));
    assert(extracted_bit.second >= 0
           && extracted_bit.first >= extracted_bit.second
           && extracted_bit.first < width);

    auto pos = dfs_walked_extract_.find(ast);
    if (pos != dfs_walked_extract_.end()
        && pos->second.find(extracted_bit) != pos->second.end()) {
      node_stack_.pop_back();
      continue;
    }
    dfs_walked_extract_[ast].insert(extracted_bit);

    smt::Op op = ast->get_op();
    if (op.is_null()) {
      if (ast->is_symbolic_const()) {
        varset_slice[ast].insert(extracted_bit);
      }
      node_stack_.pop_back();
      continue;
    } else {
      if (op.prim_op == smt::PrimOp::Ite) {
        ARG3(cond, texpr, fexpr)
        auto cond_val = solver_->get_value(cond);
        assert(cond_val->is_value());
        if (is_all_one(cond_val->to_string(), 1)) {
          node_stack_.push_back({ cond, { 0, 0 } });
          node_stack_.push_back({ texpr, extracted_bit });
        } else {
          node_stack_.push_back({ cond, { 0, 0 } });
          node_stack_.push_back({ fexpr, extracted_bit });
        }
      } else if (op.prim_op == smt::PrimOp::Implies) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        assert(unified_width(ast) == 1);
        if (!is_all_one(cond_left->to_string(), 1)) {
          node_stack_.push_back({ left, { 0, 0 } });
        } else if (is_all_one(cond_right->to_string(), 1)) {
          node_stack_.push_back({ right, { 0, 0 } });
        } else {
          node_stack_.push_back({ left, { 0, 0 } });
          node_stack_.push_back({ right, { 0, 0 } });
        }
      } else if (op.prim_op == smt::PrimOp::And) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        assert(unified_width(ast) == 1);
        if (!is_all_one(cond_left->to_string(), 1)) {
          node_stack_.push_back({ left, { 0, 0 } });
        } else if (!is_all_one(cond_right->to_string(), 1)) {
          node_stack_.push_back({ right, { 0, 0 } });
        } else {
          node_stack_.push_back({ left, { 0, 0 } });
          node_stack_.push_back({ right, { 0, 0 } });
        }
      } else if (op.prim_op == smt::PrimOp::Or) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        assert(unified_width(ast) == 1);
        if (is_all_one(cond_left->to_string(), 1)) {
          node_stack_.push_back({ left, { 0, 0 } });
        } else if (is_all_one(cond_right->to_string(), 1)) {
          node_stack_.push_back({ right, { 0, 0 } });
        } else {
          node_stack_.push_back({ left, { 0, 0 } });
          node_stack_.push_back({ right, { 0, 0 } });
        }
      } else if (op.prim_op == smt::PrimOp::BVAnd
                 || op.prim_op == smt::PrimOp::BVNand) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        std::string left_val = to_unified_bvconst(cond_left->to_string());
        std::string right_val = to_unified_bvconst(cond_right->to_string());

        // For AND, a zero bit on either side is enough; only one/one bits
        // require both operands to be tracked.
        std::vector<std::pair<int, int>> left_zero, left_one;
        find_consecutive_zeros_ones(
            left_val, { extracted_bit }, left_zero, left_one);
        for (const auto & z : left_zero) node_stack_.push_back({ left, z });

        std::vector<std::pair<int, int>> right_zero, right_one;
        find_consecutive_zeros_ones(
            right_val, left_one, right_zero, right_one);
        for (const auto & z : right_zero) node_stack_.push_back({ right, z });
        for (const auto & both_one : right_one) {
          node_stack_.push_back({ left, both_one });
          node_stack_.push_back({ right, both_one });
        }
      } else if (op.prim_op == smt::PrimOp::BVOr
                 || op.prim_op == smt::PrimOp::BVNor) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        std::string left_val = to_unified_bvconst(cond_left->to_string());
        std::string right_val = to_unified_bvconst(cond_right->to_string());

        // Dual of AND: a one bit on either side is enough; zero/zero bits
        // require both operands to be tracked.
        std::vector<std::pair<int, int>> left_zero, left_one;
        find_consecutive_zeros_ones(
            left_val, { extracted_bit }, left_zero, left_one);
        for (const auto & one : left_one) {
          node_stack_.push_back({ left, one });
        }

        std::vector<std::pair<int, int>> right_zero, right_one;
        find_consecutive_zeros_ones(
            right_val, left_zero, right_zero, right_one);
        for (const auto & one : right_one) {
          node_stack_.push_back({ right, one });
        }
        for (const auto & both_zero : right_zero) {
          node_stack_.push_back({ left, both_zero });
          node_stack_.push_back({ right, both_zero });
        }
      } else if (op.prim_op == smt::PrimOp::BVMul) {
        ARG2(left, right)
        auto cond_left = solver_->get_value(left);
        auto cond_right = solver_->get_value(right);
        assert(cond_left->is_value() && cond_right->is_value());
        std::string left_val = cond_left->to_string();
        std::string right_val = cond_right->to_string();
        auto left_w = static_cast<int>(unified_width(left));
        auto right_w = static_cast<int>(unified_width(right));
        if (is_all_zero(left_val)) {
          node_stack_.push_back({ left, { left_w - 1, 0 } });
        } else if (is_all_zero(right_val)) {
          node_stack_.push_back({ right, { right_w - 1, 0 } });
        } else {
          node_stack_.push_back({ left, { left_w - 1, 0 } });
          node_stack_.push_back({ right, { right_w - 1, 0 } });
        }
      } else if (op.prim_op == smt::PrimOp::BVAdd) {
        ARG2(left, right)
        auto msb = extracted_bit.first;
        // Carries from lower bits can affect any requested output bit.
        node_stack_.push_back({ left, { msb, 0 } });
        node_stack_.push_back({ right, { msb, 0 } });
      } else if (op.prim_op == smt::PrimOp::BVNot) {
        ARG1(back)
        node_stack_.push_back({ back, extracted_bit });
      } else if (op.prim_op == smt::PrimOp::BVXor
                 || op.prim_op == smt::PrimOp::BVXnor) {
        ARG2(left, right)
        node_stack_.push_back({ left, extracted_bit });
        node_stack_.push_back({ right, extracted_bit });
      } else if (op.prim_op == smt::PrimOp::Extract) {
        auto lsb_shift = op.idx1;
        auto msb = extracted_bit.first;
        auto lsb = extracted_bit.second;
        for (const auto & arg : *ast) {
          node_stack_.push_back(
              { arg, { static_cast<int>(lsb_shift + msb),
                       static_cast<int>(lsb_shift + lsb) } });
        }
      } else if (op.prim_op == smt::PrimOp::Concat) {
        ARG2(left, right)
        auto width_right = static_cast<int>(unified_width(right));
        auto msb = extracted_bit.first;
        auto lsb = extracted_bit.second;
        if (lsb >= width_right) {
          node_stack_.push_back(
              { left, { msb - width_right, lsb - width_right } });
        } else if (msb < width_right) {
          node_stack_.push_back({ right, { msb, lsb } });
        } else {
          node_stack_.push_back({ left, { msb - width_right, 0 } });
          node_stack_.push_back({ right, { width_right - 1, lsb } });
        }
      } else if (op.prim_op == smt::PrimOp::Zero_Extend
                 || op.prim_op == smt::PrimOp::Sign_Extend) {
        ARG1(back)
        auto width_back = static_cast<int>(unified_width(back));
        auto msb = extracted_bit.first;
        auto lsb = extracted_bit.second;
        if (lsb < width_back) {
          if (msb >= width_back) {
            node_stack_.push_back({ back, { width_back - 1, lsb } });
          } else {
            node_stack_.push_back({ back, { msb, lsb } });
          }
        }
      } else if (op.prim_op == smt::PrimOp::Equal
                 || op.prim_op == smt::PrimOp::Distinct
                 || op.prim_op == smt::PrimOp::BVComp) {
        ARG2(left, right)
        auto result_val = solver_->get_value(ast)->to_string();
        if (((op.prim_op == smt::PrimOp::Equal
              || op.prim_op == smt::PrimOp::BVComp)
             && is_all_one(result_val, 1))
            || (op.prim_op == smt::PrimOp::Distinct
                && is_all_zero(result_val))) {
          auto left_width = static_cast<int>(unified_width(left));
          auto right_width = static_cast<int>(unified_width(right));
          node_stack_.push_back({ left, { left_width - 1, 0 } });
          node_stack_.push_back({ right, { right_width - 1, 0 } });
        } else {
          // For a disequality/equality that evaluates false, one differing bit
          // is enough to justify the Boolean result.
          auto left_value = solver_->get_value(left)->to_string();
          auto right_value = solver_->get_value(right)->to_string();
          auto left_width = unified_width(left);
          auto right_width = unified_width(right);
          assert(left_width == right_width);
          assert(left_value != right_value);

          int rightmost_difference = find_rightmost_index_with_different_bit(
              left_value, right_value, left_width);
          assert(rightmost_difference != -1);

          node_stack_.push_back(
              { left, { rightmost_difference, rightmost_difference } });
          node_stack_.push_back(
              { right, { rightmost_difference, rightmost_difference } });
        }
      } else if (op.prim_op == smt::PrimOp::BVUlt
                 || op.prim_op == smt::PrimOp::BVUle
                 || op.prim_op == smt::PrimOp::BVUgt
                 || op.prim_op == smt::PrimOp::BVUge) {
        ARG2(left, right)
        auto left_value = solver_->get_value(left)->to_string();
        auto right_value = solver_->get_value(right)->to_string();
        auto left_width = unified_width(left);
        auto right_width = unified_width(right);
        assert(left_width == right_width);

        // Unsigned comparisons are decided by the most significant differing
        // bit; equal operands fall back to bit 0.
        int leftmost_difference = find_leftmost_index_with_different_bit(
            left_value, right_value, left_width);
        node_stack_.push_back(
            { left,
              { static_cast<int>(left_width) - 1, leftmost_difference } });
        node_stack_.push_back(
            { right,
              { static_cast<int>(right_width) - 1, leftmost_difference } });
      } else {
        for (const auto & arg : *ast) {
          auto arg_width = static_cast<int>(unified_width(arg));
          node_stack_.push_back({ arg, { arg_width - 1, 0 } });
        }
      }
    }
  }
}

}  // namespace pono
