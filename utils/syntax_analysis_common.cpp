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


#include "syntax_analysis_common.cpp"

#include "utils/str_util.h"
#include "utils/container_shortcut.h"

#include <cassert>

namespace pono {

namespace syntax_analysis {

// --------------------  eval_val ----------------

eval_val::eval_val(const std::string & val) {
  assert(val.find("#b") == 0);
  size_t pos = 2;
  for(; pos < val.length() ; ++ pos) {
    if ( val.at(pos) != '0' )
      break;
  }
  if (pos == val.length()) {
    // result 0
    type = type_t::NUM;
    nv = 0;
  } else {
    try {
      nv = ::pono::syntax_analysis::StrToULongLong(val.substr(pos), 2);
      type = type_t::NUM;      
    } catch (...) {
      type = type_t::STR;
      sv = val.substr(pos);
    }
  }
} // eval_val::eval_val

bool eval_val::operator<(const eval_val &r) const {
  if (type == type_t::NUM && r.type == type_t::STR)
    return true;
  if (type == type_t::STR && r.type == type_t::NUM)
    return false;
  if (type == type_t::NUM)
    return nv < r.nv;
  // both str
  if (sv.length() < r.sv.length())
    return true;
  if (sv.length() > r.sv.length())
    return false;
  for(size_t pos = 0; pos < sv.length(); ++ pos) {
    if (sv.at(pos) == '0' && r.sv.at(pos) == '1')
      return true;
    if (sv.at(pos) == '1' && r.sv.at(pos) == '0')
      return false;
  }
  return false; // equal both string, same length and save val
} // eval_val::operator<


std::string eval_val::to_string() const {
  if (type == type_t::NUM)
    return "(bv" + std::to_string(nv)+" X)";
  return "#b"+sv;
} // to_string

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

  auto term_string = new_term->to_raw_string();
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

} // namespace syntax_analysis

} // namespace pono


