/*********************                                                        */
/*! \file syntax_analysis_walker.h
 ** \verbatim
 ** Top contributors (to current version):
 **   Hongce Zhang
 ** This file is part of the pono project.
 ** Copyright (c) 2020 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief various walkers in syntax analysis
 **
 ** 
 **/

#pragma once

#include "smt-switch/smt.h"
#include "utils/syntax_analysis_common.h"

#include <set>

namespace pono {
namespace syntax_analysis {

class Walker {
protected:
  // if you want to buffer and avoid further walk
  // Skip is your chance to do that, return true
  // if already buffered
  virtual bool Skip(const smt::Term & ast) = 0;
  virtual void PreChild(const smt::Term & ast) = 0;
  virtual void PostChild(const smt::Term & ast) = 0;
public:
  void WalkDFS(const smt::Term & ast);
  void WalkBFS(const smt::Term & ast);
}; // Walker


// -----------------------------------------------------

// we'd better extract from msat's term
// btor's mixing bool and 1-width bv will
// create problems
class OpExtractor : public Walker {

public:
  syntax_analysis::SyntaxStructure & 
    GetSyntaxConstruct() {
      return constructs; }

protected:
  std::unordered_set<smt::Term> walked_nodes_;
  std::unordered_set<smt::Term> all_symbols_;
  syntax_analysis::SyntaxStructure constructs;
  
  virtual bool Skip(const smt::Term & ast) override;
  virtual void PreChild(const smt::Term & ast) override;
  virtual void PostChild(const smt::Term & ast) override;


}; // class OpExtractor

// -----------------------------------------------------

class TermExtractor : public Walker {
public:
  // ----------- TYPE --------------- //
  struct node_info_t {
    bool in;
    bool related;
    unsigned level;
    node_info_t() : in(false), related(false), level(0) {}
  };
  
  // ----------- CONSTRUCTOR --------------- //
  // if level > 0, then we will count level, otherwise, we don't care about the levels
  TermExtractor(const std::unordered_set<smt::Term> & varset, bool collect_constants, unsigned level,
      std::map<unsigned, std::vector<smt::Term>> & width_to_term_table,
      std::unordered_set<smt::Term> & all_terms_set
      //std::unordered_set<smt::Term> & related_terms_set
      ) :
    related_vars_(varset), collect_constants_(collect_constants), level_(level),
    width_to_terms_(width_to_term_table), all_terms_(all_terms_set)
    //related_terms_(related_terms_set) 
    { }
    
  //const std::map<unsigned, std::vector<smt::Term>> & GetTermsByWidth() const { return width_to_terms_; }
  //const std::unordered_set<smt::Term> & GetAllTerms() const { return all_terms_; }
  const std::map<unsigned, std::vector<smt::Term>> & GetConstants() const { return width_to_constants_; }
  
  // public method inherited: WalkDFS (*) /BFS
  
protected:
  std::unordered_map<smt::Term, node_info_t> walked_nodes_;
  const std::unordered_set<smt::Term> & related_vars_; // we will also bring in unrelated vars
  bool collect_constants_;
  unsigned level_;
  
  std::map<unsigned, std::vector<smt::Term>> & width_to_terms_;
  std::map<unsigned, std::vector<smt::Term>> width_to_constants_; // const is not needed, as you may not always need it
  std::unordered_set<smt::Term> & all_terms_;
  // std::unordered_set<smt::Term> & related_terms_;
  
  virtual bool Skip(const smt::Term & ast) override;
  virtual void PreChild(const smt::Term & ast) override;
  virtual void PostChild(const smt::Term & ast) override;

}; // class TermExtractor

// -----------------------------------------------------


class TermScore : public Walker {
public:
  // ----------- TYPE --------------- //
  struct term_score_t {
    unsigned score;
    term_score_t(unsigned s) : score(s) {}
  };
  
  typedef std::unordered_map<smt::Term, term_score_t> score_map_t;

  TermScore() {} // do nothing
  const score_map_t & GetScoreMap() const {return scores_;}
  
protected:
  score_map_t scores_;
  
  virtual bool Skip(const smt::Term & ast) override;
  virtual void PreChild(const smt::Term & ast) override;
  virtual void PostChild(const smt::Term & ast) override;

}; // TermScore


// -----------------------------------------------------

class ParentExtract : public Walker {
public:
  // ----------- TYPE --------------- //
  typedef std::unordered_map<smt::Term, smt::UnorderedTermSet> parent_map_t;

  ParentExtract() {} // do nothing
  void ClearCache() { parent_.clear(); }
  const parent_map_t & GetParentRelation() {return parent_;}
  bool RegisterNewParentRelation(const smt::Term &child, const smt::Term &parent) {
    auto ret = parent_[child].insert(parent);
    return ret.second;
  }
  
protected:

  std::unordered_set<smt::Term> walked_nodes_;
  parent_map_t parent_;
  
  virtual bool Skip(const smt::Term & ast) override;
  virtual void PreChild(const smt::Term & ast) override;
  virtual void PostChild(const smt::Term & ast) override;

}; // ParentExtract

// -----------------------------------------------------

class ConstantExtractor: public Walker {
public:
  // ----------- TYPE --------------- //
  typedef std::map<unsigned, std::vector<smt::Term>> width_constant_map_t;

  ConstantExtractor(width_constant_map_t & out,
    std::unordered_set<std::string> & cnstr_strs
    ) : width_constant_map(out), constants_strs_(cnstr_strs)  {}

protected:
  width_constant_map_t & width_constant_map;
  std::unordered_set<std::string> & constants_strs_;
  std::unordered_set<smt::Term> walked_nodes_;

  virtual bool Skip(const smt::Term & ast) override;
  virtual void PreChild(const smt::Term & ast) override;
  virtual void PostChild(const smt::Term & ast) override;

}; // ConstantExtractor

// -----------------------------------------------------

class SliceExtractor: public Walker {
public:
  using width_term_map_t = PerVarsetInfo::width_term_map_t;
  //typedef std::map<unsigned, PerWidthInfo> width_terms_map_t;
  typedef std::pair<unsigned, unsigned> ext_position_t;
  typedef std::unordered_map<smt::Term, std::set<ext_position_t>> sv2exts_t;

  SliceExtractor(
    /*OUTPUT*/ width_term_map_t & ext_terms,
    /*INPUT*/  const std::unordered_set<smt::Term> & varset) :
    ext_terms_(ext_terms), related_vars_(varset) { }
  const sv2exts_t & GetSvSlice() const { return sv2exts_; }

protected:
  std::unordered_set<smt::Term> walked_nodes_;
  width_term_map_t & ext_terms_; // point to width 1
  sv2exts_t sv2exts_;
  const std::unordered_set<smt::Term> & related_vars_; 

  virtual bool Skip(const smt::Term & ast) override;
  virtual void PreChild(const smt::Term & ast) override;
  virtual void PostChild(const smt::Term & ast) override;

}; // SliceExtractor

// -----------------------------------------------------

// you may also want to register the model -> full model map
class TermLearner {

public:
  TermLearner( to_next_t to_next_func, 
      score_t score_func,
      //cex_term_map_t & cex_pred_map, 
      smt::SmtSolver & btor, ParentExtract & parent_extractor) : 
    to_next_(to_next_func),
    score_(score_func),
    // cex_pred_map_ref_(cex_pred_map),
    solver_(btor),
    parent_extractor_(parent_extractor) {}
    
  unsigned learn_terms_from_cex(IC3FormulaModel * pre, IC3FormulaModel * post, 
    const smt::Term & trans,
    /*OUTPUT*/  PerVarsetInfo & varset_info );

  unsigned learn_terms_from_cex_same_frame(IC3FormulaModel * pre, IC3FormulaModel * post, 
    const smt::Term & trans,
    /*OUTPUT*/  PerVarsetInfo & varset_info );

  unsigned vars_extract_bit_level(IC3FormulaModel * post,  /*OUTPUT*/  PerVarsetInfo & varset_info) ;
  
protected:
  // smt::Term trans_; // you need to give the full constraints and etc...
  to_next_t to_next_;
  score_t score_;
  //cex_term_map_t & cex_pred_map_ref_; // from unsat enum
  smt::SmtSolver & solver_;
  ParentExtract & parent_extractor_; // ParentExtract

  unsigned syntactic_score_factor = 2;
  unsigned syntactic_score_delta = 0;
  
protected:
  unsigned same_val_replace_ast( /*INOUT*/  PerVarsetInfo & varset_info );
  unsigned same_val_replace_ast_same_frame( /*INOUT*/  PerVarsetInfo & varset_info );
  unsigned replace_hierachically(
    const smt::Term & orig, const smt::Term & repl, /*INOUT*/  PerVarsetInfo & varset_info );
  
  unsigned replace_hierachically_w_parent(
    const smt::Term & orig, const smt::Term & repl, PerVarsetInfo & varset_info,
    smt::TermVec & output_new_terms ) ;

  unsigned concat_to_extract(/*INOUT*/  PerVarsetInfo & varset_info);
  unsigned extract_complement(/*INOUT*/  PerVarsetInfo & varset_info);

}; // class TermLearner

} // namespace syntax_analysis
} // namespace pono

