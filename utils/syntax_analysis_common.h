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
#include <functional>

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
struct eval_val {
  // the starting char is always 1
  std::string sv;
  
  eval_val(const std::string & val); 
  // will remove #b0...0 and then decide to convert or not
  // default copy and assignment, and then

  bool operator==(const eval_val &r) const {
    return (sv == r.sv) ;
  }

  bool operator<(const eval_val &) const;

  std::string to_string() const {return sv;}

  // the first one is always 1....
  // so, if one is shorter, it must be smaller

}; // struct eval_val


struct eval_val_hash {
  std::size_t operator() (const eval_val & k) const {
    return (std::hash<std::string>()(k.sv));
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
  std::unordered_set<std::string> predicates_str;
  // std::unordered_map<smt::Term, smt::Term> pred_next_to_pred_curr;
  const PerVarsetInfo & varset_info; // reference from VarTermManager
  unsigned prev_refine_constraint_count;
  std::map<unsigned, term_const_num> prev_per_width_term_num;

  PerCexInfo(const PerVarsetInfo & info,
    unsigned constraint_count) : varset_info(info), prev_refine_constraint_count(constraint_count) {}
  void ResetPredicates() { predicates_nxt.clear(); predicates_str.clear(); }
};

// cex -> related information
// Here I don't use `IC3Formula *` because it is changing
// and not reliable
typedef std::unordered_map<IC3FormulaModel *, PerCexInfo>   cex_term_map_t; // the enumeration position of a cex

// --------------------------------------------------
//
//               Syntax Structures
//
// --------------------------------------------------

  std::string name_sanitize(const std::string &); 
  std::string name_desanitize(const std::string &s);
  uint64_t get_width_of_var(const smt::Term & v);
  smt::Term smt_string_to_const_term(const std::string & val, smt::SmtSolver & s);
  bool is_primop_symmetry(smt::PrimOp);

  struct sygus_op {
    virtual smt::Op to_smt_op() const = 0;
  };

  struct concat_t : public sygus_op {
    uint64_t width1; uint64_t width2;
    concat_t(uint64_t w1, uint64_t w2) :
      width1(w1), width2(w2) { }
    
    virtual smt::Op to_smt_op() const override {
      return smt::Op(smt::PrimOp::Concat);
    }

  };
  struct extract_t : public sygus_op {
    uint64_t input_width, h, l;
    extract_t(uint64_t iw, uint64_t _h, uint64_t _l) :
      input_width(iw), h(_h), l(_l) { }
    virtual smt::Op to_smt_op() const override {
      return smt::Op(smt::PrimOp::Extract, h, l);
    }
  };
  struct rotate_t : public sygus_op {
    smt::PrimOp op; uint64_t param;
    rotate_t(smt::PrimOp _op, uint64_t _param) :
      op(_op), param(_param) { }
    virtual smt::Op to_smt_op() const override {
      return smt::Op(op, param);
    }
  };
  struct extend_t : public sygus_op {
    smt::PrimOp op; uint64_t extw, input_width;
    extend_t(smt::PrimOp _op, uint64_t _extw, uint64_t _iw) :
      op(_op), extw(_extw), input_width(_iw) { }
    virtual smt::Op to_smt_op() const override {
      return smt::Op(op, extw);
    }
  };

  class concat_hash {
  public:
    size_t operator() (const concat_t & t) const;
  }; 
  class extract_hash {
  public:
    size_t operator() (const extract_t & t) const;
  }; 
  class rotate_hash {
  public:
    size_t operator() (const rotate_t & t) const;
  }; 
  class extend_hash {
  public:
    size_t operator() (const extend_t & t) const;
  }; 


  bool operator==(const concat_t & a,  const concat_t & b);
  bool operator==(const extract_t & a, const extract_t & b);
  bool operator==(const rotate_t & a,  const rotate_t & b);
  bool operator==(const extend_t & a,  const extend_t & b);

  struct BvConstructs {
    std::unordered_set<std::string> symbol_names;
    // let's use to_string to fill it? so we hope we don't need to add | ourselves
    std::unordered_set<std::string> constants; // let's convert it to string
    std::unordered_set<smt::PrimOp> op_unary;  // unary operators: (_ bv x) -> (_ bv x)
    std::unordered_set<smt::PrimOp> op_binary; // binary operators: (_ bv x) x (_ bv x) -> (_ bv x)
    std::unordered_set<smt::PrimOp> op_comp;  // comparators: (_ bv x) x (_ bv x) -> bool
    // set of (width1, width2)
    std::unordered_set<concat_t, concat_hash> op_concat;
    // set of (input_width, h, l)
    std::unordered_set<extract_t, extract_hash> op_extract;
    // set of (op, param)
    std::unordered_set<rotate_t, rotate_hash> op_rotate;
    // set of (op, param, input_width)
    std::unordered_set<extend_t, extend_hash> op_extend;

    // default constructor
    BvConstructs() {}

  }; // class BvConstructs

  class SyntaxStructure{

  public:
    typedef std::unordered_map<uint64_t, BvConstructs> SyntaxT;

    bool new_constructs;

    const SyntaxT & GetSyntaxConstruct() const {
        return syntax_; }
    
    void insert_symbol   (uint64_t width, const std::string & name);
    void insert_const    (uint64_t width, const std::string & name);
    void insert_op_unary (uint64_t width, smt::PrimOp);
    void insert_op_binary(uint64_t width, smt::PrimOp);
    void insert_op_comp  (uint64_t width, smt::PrimOp);
    void insert_concat   (uint64_t width, concat_t  && );
    void insert_extract  (uint64_t width, extract_t && );
    void insert_rotate   (uint64_t width, rotate_t  && );
    void insert_extend   (uint64_t width, extend_t  && );

    void CutVars(
      const std::unordered_set<std::string> & keep_vars_name,
      const std::unordered_set<std::string> & remove_vars_name);

    void RemoveUnusedStructure();
    void RemoveExtract();
    void RemoveConcat();
    void AddBvultBvule();
    void AndOrConvert();

    static smt::Term const_to_term(const std::string & val, smt::SmtSolver & s);
    
  protected:
    SyntaxT syntax_;

  }; // 

  typedef std::unordered_map<uint64_t, BvConstructs> SyntaxStructureT;

  typedef std::function<smt::Term(const smt::Term &)> to_next_t; 
  typedef std::function<unsigned(const smt::Term &)> score_t; 
  // function signature that can convert convert ast to next

} // namespace syntax_analysis
} // namespace pono


