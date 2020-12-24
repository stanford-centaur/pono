/*********************                                                        */
/*! \file syntax_analysis_common.h
 ** \verbatim
 ** Top contributors (to current version):
 **   Hongce Zhang
 ** This file is part of the pono project.
 ** Copyright (c) 2020 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief Various common data structures for syntax analysis
 **
 ** 
 **/
 
 
#pragma once

#include "smt-switch/smt.h"
#include "utils/sygus_ic3formula_helper.h"

#include <map>

namespace pono {

namespace syntax_analysis {


// record what terms and constants are available
// in the original problem
struct PerWidthInfo {
 std::vector<smt::Term> terms;
 std::vector<smt::Term> constants;
};

class VarTermManager;


struct PerVarsetInfo {
  // --- type definition --- //
  typedef std::map<unsigned, smt::TermVec> width_term_map_t;
  struct state_t  {
    enum stage_t {EMPTY, WPARTIAL, WALL, FROMCEX, EXTRACTBITS, VCLT, VCLTE} stage;
    unsigned partial_width_done;
    // you don't need to cache those constants, already in width_to_constants_
    //state_t() : stage(EMPTY), partial_width_done(0) {}
    explicit state_t(stage_t s) : stage(s), partial_width_done(0)  {}
  }; // class state_t

  std::map<unsigned, PerWidthInfo> terms;
  std::unordered_set<std::string> terms_strings;
  // --- more info --- //

  //PerVarsetInfo() {}
  PerVarsetInfo(state_t::stage_t s): state(s) {}

  // if true, inserting is done
  bool TermLearnerInsertTerm(const smt::Term & t);
  bool TermLearnerIsOut(const smt::Term & p) const {
    return all_terms.find(p) == all_terms.end(); }

  bool use_lt() const { return state.stage == state_t::VCLT; }
  bool use_lte() const { return state.stage == state_t::VCLTE; }

protected: // let's restrict the accesses to these fields
  friend class VarTermManager; 
  state_t state;
  width_term_map_t terms_buffer; // will be copied from term_walker
  smt::UnorderedTermSet all_terms;
}; // class PerVarsetInfo


// --------------------------------------------------
//
//               eval_val
//
// --------------------------------------------------

// value for enumeration
struct eval_val { // will always convert to uint64_t, if width < 64

  enum type_t {NUM, STR} type;
  uint64_t nv;
  std::string sv;
  
  eval_val(const std::string & val); // will remove #b0...0 and then decide to convert or not
  // default copy and assignment, and then

  bool operator==(const eval_val &r) const {
    return (type == r.type) && 
      (type != type_t::NUM || nv == r.nv) && 
      (type != type_t::STR || sv == r.sv);
  }

  bool operator<(const eval_val &) const;

  std::string to_string() const;

  // the first one is always 1....
  // so, if one is shorter, it must be smaller

}; // struct eval_val


struct eval_val_hash {
  std::size_t operator() (const eval_val & k) const {
    return (k.type == k.NUM ? std::hash<uint64_t>()(k.nv) : std::hash<std::string>()(k.sv));
  }
};

// --------------------------------------------------
//
//               PerCexInfo
//
// --------------------------------------------------

struct PerCexInfo {
  struct term_const_num{
    unsigned term_num;
    unsigned const_num;
    term_const_num(): term_num(0), const_num(0) {}
  };

  std::unordered_map<smt::Term,eval_val> terms_val_under_cex;
  std::vector<smt::Term> predicates_nxt;
  std::unordered_map<smt::Term, smt::Term> pred_next_to_pred_curr;
  const PerVarsetInfo & varset_info; // reference from VarTermManager

  std::map<unsigned, term_const_num> prev_per_width_term_num;
  PerCexInfo(const PerVarsetInfo & info) : varset_info(info) {}
};

// cex -> related information
// Here I don't use `IC3Formula *` because it is changing
// and not reliable
typedef std::unordered_map<IC3FormulaModel *, PerCexInfo>   cex_term_map_t; // the enumeration position of a cex


} // namespace syntax_analysis
} // namespace pono


