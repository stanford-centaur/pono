/*********************                                                        */
/*! \file syntax_analysis_common.cpp
 ** \verbatim
 ** Top contributors (to current version):
 **   Hongce Zhang
 ** This file is part of the pono project.
 ** Copyright (c) 2020 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief Some type defs, shared by apdr/enum/learner
 **
 ** 
 **/


#include "utils/syntax_analysis_common.h"
#include "utils/str_util.h"
#include "utils/container_shortcut.h"

#include <cassert>
#include <queue>

namespace pono {

namespace syntax_analysis {



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
  return true;
}


std::string convert_bin_str_to_decimal(const std::string & in) {
  std::vector<char> out = {0};
  for (auto pos = in.begin(); pos != in.end(); ++pos) {
    syntax_analysis::mul2(out); // initially it is 0 so okay
    if (*pos == '1')
      syntax_analysis::add1(out);
  }
  std::string ret;
  for (auto pos = out.rbegin(); pos != out.rend(); ++pos) {
    ret += *pos + '0';
  }
  return ret;
}


// --------------------  eval_val ----------------

eval_val::eval_val(const std::string & val) {
  if (val == "true" || val == "false") {
    sv= (val == "true") ? "1" : "0";
    return;
  }

  if(val.find("#b") != 0) { // then it is (_ bvX width)
    std::string width; // width is no use
    bool succ = extract_decimal_width(val, sv, width);
    assert(succ);
    return;
  }
  
  // #b something here
  size_t pos = 2;
  for(; pos < val.length() ; ++ pos) {
    if ( val.at(pos) != '0' )
      break;
  }
  if (pos == val.length())  {
    sv = "0";
    return;
  }
  // convert 101 --> 5
  sv = convert_bin_str_to_decimal(val.substr(pos));
} // eval_val::eval_val



bool eval_val::operator<(const eval_val &r) const {
  // both str
  if (sv.length() < r.sv.length())
    return true;
  if (sv.length() > r.sv.length())
    return false;
  for(size_t pos = 0; pos < sv.length(); ++ pos) {
    //if (sv.at(pos) == '0' && r.sv.at(pos) == '1')
    if (sv.at(pos) < r.sv.at(pos))
      return true;
    //if (sv.at(pos) == '1' && r.sv.at(pos) == '0')
    if (sv.at(pos) > r.sv.at(pos))
      return false;
  }
  return false; // equal both string, same length and save val
} // eval_val::operator<


// --------------------  PerVarsetInfo ----------------
static unsigned get_width(const smt::Term & t) {
  auto sort_kind = t->get_sort()->get_sort_kind() ;
  if ( sort_kind == smt::SortKind::BOOL)
    return 1; // also make it bv?
  else if (sort_kind == smt::SortKind::BV)
    return t->get_sort()->get_width();
  assert(false); // we don't know how to handle
  return 0;
}

// if true, inserting is done
bool PerVarsetInfo::TermLearnerInsertTerm(const smt::Term & new_term) {

  auto term_string = new_term->to_string();
  auto ins_res = terms_strings.insert(term_string);
  if (!ins_res.second) // if already exists, will not insert
    return false;

  unsigned width = get_width(new_term);
  assert (IN(width, terms)); // concat -> extract does not change this
  bool is_val = (new_term->is_value());
  if(is_val)
    terms.at(width).constants.push_back(new_term);
  else
    terms.at(width).terms.push_back(new_term);
  all_terms.insert(new_term);

  return true;
} // PerVarsetInfo::TermLearnerInsertTerm


// --------------------------------------------------
//
//               Syntax Structures
//
// --------------------------------------------------


std::string name_sanitize(const std::string &s) {
  if (s.length() > 2 && s.front() == '|' && s.back() == '|')
    return s; // already | |
  bool need_separator = false;
  for (auto pos = s.begin(); pos != s.end(); ++pos) {
    char c = *pos;
    if (isalnum(c))
      continue;
    if (c == '.' || c == '-' )
      continue;
    if (c == '_' && pos != s.begin())
      continue;
    // else
    need_separator = true;
    break;
  }
  if (need_separator)
    return "|" + s + "|";
  return s;
}

std::string name_desanitize(const std::string &s) {
  if (s.length() > 2 && s.front() == '|' && s.back() == '|')
    return s.substr(1,s.length()-2); // already | |
  return s;
}

uint64_t get_width_of_var(const smt::Term & v) {
  if ( v->get_sort()->get_sort_kind() == smt::SortKind::BOOL )
    return 0;
  else if ( v->get_sort()->get_sort_kind() == smt::SortKind::BV)
   return v->get_sort()->get_width();
  assert (false);
  return 0;
}


smt::Term smt_string_to_const_term(const std::string & val, smt::SmtSolver & s) {
  if (val == "true" || val == "True") {
    return s->make_term(true);
  } else if (val == "false" || val == "False") {
    return s->make_term(false);
  } else if(val.find("#b") == 0) {
    uint64_t width = val.length()-2; // exclude #b
    auto bvsort = s->make_sort(smt::SortKind::BV, width);
    return s->make_term(val.substr(2), bvsort, 2UL);
  } else if(val.find("#x") == 0) {
    uint64_t width = (val.length()-2)*4; // exclude #x
    auto bvsort = s->make_sort(smt::SortKind::BV, width);
    return s->make_term(val.substr(2), bvsort, 16UL);
  } else if (val.find("(_ bv") == 0) {
    auto v_w_pair = SplitSpaceTabEnter(val.substr(5));
    assert(v_w_pair.size()>=2);
    uint64_t width = (uint64_t)(StrToULongLong(ReplaceAll(v_w_pair[1],")",""), 10));
    auto bvsort = s->make_sort(smt::SortKind::BV, width);
    return s->make_term(v_w_pair[0], bvsort, 10UL);    
  } 
  assert (false); // unknown format
  return nullptr;
}


bool is_primop_symmetry(smt::PrimOp op) {
    switch (op) {
      case smt::PrimOp::And:
      case smt::PrimOp::Or:
      case smt::PrimOp::Xor:
      // case smt::PrimOp::Iff: why?
      case smt::PrimOp::BVAnd:
      case smt::PrimOp::BVOr:
      case smt::PrimOp::BVXor:
      case smt::PrimOp::BVNor:
      case smt::PrimOp::BVXnor:
      case smt::PrimOp::BVNand:
      case smt::PrimOp::BVAdd:
      case smt::PrimOp::BVMul:
      case smt::PrimOp::Equal:
      case smt::PrimOp::Distinct: return true;
/*
      case smt::PrimOp::Not:
      case smt::PrimOp::BVNeg:
      case smt::PrimOp::BVNot:
      case smt::PrimOp::Implies: 
      case smt::PrimOp::BVSub:
      case smt::PrimOp::BVUdiv:
      case smt::PrimOp::BVSdiv:
      case smt::PrimOp::BVUrem:
      case smt::PrimOp::BVSrem:
      case smt::PrimOp::BVSmod:
      case smt::PrimOp::BVShl:
      case smt::PrimOp::BVAshr:
      case smt::PrimOp::BVLshr: 
      case smt::PrimOp::BVComp:
      case smt::PrimOp::BVUlt:
      case smt::PrimOp::BVUle:
      case smt::PrimOp::BVUgt:
      case smt::PrimOp::BVUge:
      case smt::PrimOp::BVSlt:
      case smt::PrimOp::BVSle:
      case smt::PrimOp::BVSgt:
      case smt::PrimOp::BVSge: 
      case .. others*/
      default:
        return false;
  } // switch
  assert (false);
} // is_primop_symmetry

// ------------------------------------------------------

bool operator==(const concat_t & a, const concat_t & b) {
  return ( (a.width1 == b.width1) && (a.width2 == b.width2) );
}
bool operator==(const extract_t & a, const extract_t & b) {
  return ( (a.h == b.h) && (a.l == b.l) && (a.input_width == b.input_width) );
}
bool operator==(const rotate_t & a, const rotate_t & b) {
  return ( (a.op == b.op) && (a.param == b.param) );
}
bool operator==(const extend_t & a, const extend_t & b) {
  return ( (a.input_width == b.input_width) && (a.op == b.op) && (a.extw == b.extw) );
}

template <class T>
inline size_t hash_combine(std::size_t seed, const T& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    return seed;
}

size_t concat_hash::operator() (const concat_t & t) const {
  std::hash<uint64_t> hasher;
  return hash_combine(hasher(t.width2), t.width2);
}

size_t extract_hash::operator() (const extract_t & t) const {
  std::hash<uint64_t> hasher;
  return hash_combine(hash_combine(hasher(t.h), t.l), t.input_width);
}

size_t rotate_hash::operator() (const rotate_t & t) const {
  std::hash<uint64_t> hasher;
  return hash_combine(hasher(t.param), t.op);
}

size_t extend_hash::operator() (const extend_t & t) const {
  std::hash<uint64_t> hasher;
  return hash_combine(hash_combine(hasher(t.input_width), t.extw), t.op);
}

// ------------------------------------------------------

    
void SyntaxStructure::insert_symbol   (uint64_t width, const std::string & name) {
  if ( !IN(width, syntax_) ) {
    syntax_.emplace(width, BvConstructs());
    new_constructs = true;
  }

  BvConstructs & cnstr = syntax_.at(width);
  auto res = cnstr.symbol_names.insert( name_sanitize ( name ));
  if(res.second)
    new_constructs = true; // indeed inserted

}

void SyntaxStructure::insert_const    (uint64_t width, const std::string & val) {
  if ( !IN(width, syntax_) ) {
    syntax_.emplace(width, BvConstructs());
    new_constructs = true;
  }
  BvConstructs & cnstr = syntax_.at(width);

  auto res = cnstr.constants.insert(val);
  if(res.second)
    new_constructs = true; // indeed inserted
}

void SyntaxStructure::insert_op_unary (uint64_t width, smt::PrimOp op) {
  if ( !IN(width, syntax_) ) {
    syntax_.emplace(width, BvConstructs());
    new_constructs = true;
  }
  BvConstructs & cnstr = syntax_.at(width);
  auto res = cnstr.op_unary.insert(op); 
  new_constructs = res.second ? true: new_constructs;
}

void SyntaxStructure::insert_op_binary(uint64_t width, smt::PrimOp op) {
  if ( !IN(width, syntax_) ) {
    syntax_.emplace(width, BvConstructs());
    new_constructs = true;
  }
  BvConstructs & cnstr = syntax_.at(width);
  auto res = cnstr.op_binary.insert(op); 
  new_constructs = res.second ? true: new_constructs;
}

void SyntaxStructure::insert_op_comp  (uint64_t width, smt::PrimOp op) {
  if ( !IN(width, syntax_) ) {
    syntax_.emplace(width, BvConstructs());
    new_constructs = true;
  }
  BvConstructs & cnstr = syntax_.at(width);
  auto res = cnstr.op_comp.insert(op); 
  new_constructs = res.second ? true: new_constructs;
}

void SyntaxStructure::insert_concat   (uint64_t width, concat_t && c ) {
  if ( !IN(width, syntax_) ) {
    syntax_.emplace(width, BvConstructs());
    new_constructs = true;
  }
  BvConstructs & cnstr = syntax_.at(width);
  auto res = cnstr.op_concat.insert(c);
  new_constructs = res.second ? true: new_constructs;
}

void SyntaxStructure::insert_extract  (uint64_t width, extract_t && c) {
  if ( !IN(width, syntax_) ) {
    syntax_.emplace(width, BvConstructs());
    new_constructs = true;
  }
  BvConstructs & cnstr = syntax_.at(width);
  auto res = cnstr.op_extract.insert(c);
  new_constructs = res.second ? true: new_constructs;
}

void SyntaxStructure::insert_rotate   (uint64_t width, rotate_t && c) {
  if ( !IN(width, syntax_) ) {
    syntax_.emplace(width, BvConstructs());
    new_constructs = true;
  }
  BvConstructs & cnstr = syntax_.at(width);
  auto res = cnstr.op_rotate.insert(c);
  new_constructs = res.second ? true: new_constructs;
}

void SyntaxStructure::insert_extend   (uint64_t width, extend_t && c) {
  if ( !IN(width, syntax_) ) {
    syntax_.emplace(width, BvConstructs());
    new_constructs = true;
  }
  BvConstructs & cnstr = syntax_.at(width);
  auto res = cnstr.op_extend.insert(c);
  new_constructs = res.second ? true: new_constructs;
}


void SyntaxStructure::CutVars(
  const std::unordered_set<std::string> & keep_vars_name,
  const std::unordered_set<std::string> & remove_vars_name) {
  
  for (auto & w_cnstr : syntax_) {
    auto width = w_cnstr.first;
    auto & cnstr = w_cnstr.second;
    auto pos = cnstr.symbol_names.begin();
    while(pos != cnstr.symbol_names.end()) {
      const auto & name = *pos;
      if (!keep_vars_name.empty() && (
           !IN(name ,keep_vars_name ) &&
           !IN(name_sanitize(name) ,keep_vars_name ) && 
           !IN(name_desanitize(name) ,keep_vars_name ) ) ) {
        // cut it
        pos = cnstr.symbol_names.erase(pos);
        continue;
      }
      if (IN(name, remove_vars_name) || 
          IN(name_sanitize(name),remove_vars_name) ||
          IN(name_desanitize(name),remove_vars_name) ) { 
        pos = cnstr.symbol_names.erase(pos);
        continue;
      }
      ++ pos;
    }
  } // here we compute the vars to keep
}

void SyntaxStructure::RemoveExtract() {
  for (auto & width_cnstr : syntax_) {
    auto width = width_cnstr.first;
    auto & cnstr = width_cnstr.second;
    cnstr.op_extract.clear();
  }
}
void SyntaxStructure::RemoveConcat() {
  for (auto & width_cnstr : syntax_) {
    auto width = width_cnstr.first;
    auto & cnstr = width_cnstr.second;
    cnstr.op_concat.clear();
  }
}

void SyntaxStructure::AddBvultBvule() {
  for (auto & width_cnstr : syntax_) {
    auto width = width_cnstr.first;
    auto & cnstr = width_cnstr.second;
    if (width == 0 || width == 1) {
      // remove BVUlt BVUle 
      cnstr.op_comp.erase(smt::PrimOp::Equal);
      cnstr.op_comp.erase(smt::PrimOp::Distinct);
      cnstr.op_comp.erase(smt::PrimOp::BVUlt);
      cnstr.op_comp.erase(smt::PrimOp::BVUle);
      cnstr.op_comp.erase(smt::PrimOp::BVUgt);
      cnstr.op_comp.erase(smt::PrimOp::BVUge);
      cnstr.op_comp.erase(smt::PrimOp::BVSlt);
      cnstr.op_comp.erase(smt::PrimOp::BVSle);
      cnstr.op_comp.erase(smt::PrimOp::BVSgt);
      cnstr.op_comp.erase(smt::PrimOp::BVSge);
    } else {
      cnstr.op_comp.insert(smt::PrimOp::BVUlt);
      cnstr.op_comp.insert(smt::PrimOp::BVUle);
      cnstr.op_comp.erase(smt::PrimOp::BVUgt);
      cnstr.op_comp.erase(smt::PrimOp::BVUge);
      cnstr.op_comp.erase(smt::PrimOp::BVSlt);
      cnstr.op_comp.erase(smt::PrimOp::BVSle);
      cnstr.op_comp.erase(smt::PrimOp::BVSgt);
      cnstr.op_comp.erase(smt::PrimOp::BVSge);
      // find and remove others
    }
  }
}

void SyntaxStructure::AndOrConvert() {
  for (auto & width_cnstr : syntax_) {
    auto width = width_cnstr.first;
    auto & cnstr = width_cnstr.second;
    if ( ( IN( smt::PrimOp::BVAnd , cnstr.op_binary ) || IN( smt::PrimOp::And , cnstr.op_binary ) )&& 
         ( IN( smt::PrimOp::BVOr , cnstr.op_binary )  || IN( smt::PrimOp::Or , cnstr.op_binary ) ) && 
         ( IN( smt::PrimOp::BVNot , cnstr.op_unary ) || IN( smt::PrimOp::Not , cnstr.op_unary ) ) ) {
      cnstr.op_binary.erase(smt::PrimOp::BVOr);
      cnstr.op_binary.erase(smt::PrimOp::Or);
    }
  }
}
  
void SyntaxStructure::RemoveUnusedStructure() {
 // todo 
  typedef uint64_t width_t;
  //std::unordered_set<width_t> start_set;
  std::unordered_map<width_t, std::unordered_set<width_t>> use_map;
  std::unordered_map<width_t, std::unordered_set<width_t>> to_map;

  std::unordered_set<width_t> reachable_set;
  std::queue<width_t> q;


  for (auto && width_cnstr : syntax_) {
    auto width = width_cnstr.first;
    const auto & cnstr = width_cnstr.second;
    if (!cnstr.constants.empty() || 
        !cnstr.symbol_names.empty()) {
      q.push(width);
      use_map[width].insert(width);
    }
    if (!cnstr.op_unary.empty() || !cnstr.op_binary.empty() || !cnstr.op_concat.empty())
      use_map[width].insert(width);
    for (const auto & exd: cnstr.op_extend)
      use_map[width].insert(exd.input_width);
    for (const auto & extract: cnstr.op_extract)
      use_map[width].insert(extract.input_width);
    for (const auto & concats: cnstr.op_concat) {
      use_map[width].insert(concats.width1);      
      use_map[width].insert(concats.width2);      
    }
    // rotate is between same width
  }

  // revert use_map to dst_map
  for (auto && dst_srcs_pair : use_map) {
    for (auto && src : dst_srcs_pair.second) {
      to_map[src].insert(dst_srcs_pair.first);
    }
  }

  // now do the graph reach algorithm
  while(!q.empty()) {
    width_t cur = q.front();
    q.pop();

    for (auto dstw : to_map[cur]) {
      if(!IN(dstw,reachable_set )) {
        reachable_set.insert(dstw);
        q.push(dstw);
      }
    }
  } // while queue is not empty
  auto pos = syntax_.begin();
  while (pos != syntax_.end()) {
    if (IN(pos->first, reachable_set)) {
      // for the uses, remove exd, extract concat
      { // exd
        auto exd_pos = pos->second.op_extend.begin();
        while (exd_pos != pos->second.op_extend.end())  {
          if (!IN(exd_pos->input_width, reachable_set))
            exd_pos = pos->second.op_extend.erase(exd_pos);
          else
            ++ exd_pos;
        }
      } // end of exd
      { // extract
        auto extract_pos = pos->second.op_extract.begin();
        while (extract_pos != pos->second.op_extract.end())  {
          if (!IN(extract_pos->input_width, reachable_set))
            extract_pos = pos->second.op_extract.erase(extract_pos);
          else
            ++ extract_pos;
        }
      } // end of extract
      { // concat
        auto concat_pos = pos->second.op_concat.begin();
        while (concat_pos != pos->second.op_concat.end())  {
          if (!IN(concat_pos->width1, reachable_set) || !IN(concat_pos->width2, reachable_set))
            concat_pos = pos->second.op_concat.erase(concat_pos);
          else
            ++ concat_pos;
        }
      } // end of concat
      ++pos;
    }
    else
      pos = syntax_.erase(pos);
  } // remove not in reachable set
}


} // namespace syntax_analysis
} // namespace pono


