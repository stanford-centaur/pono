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

class ParentExtract : public Walker {
public:
  // ----------- TYPE --------------- //
  typedef std::unordered_map<smt::Term, smt::UnorderedTermSet> parent_map_t;

  ParentExtract() {} // do nothing
  static void ClearCache() { parent_.clear(); }
  static const parent_map_t & GetParentRelation() {return parent_;}
  static bool RegisterNewParentRelation(const smt::Term &child, const smt::Term &parent) {
    auto ret = parent_[child].insert(parent);
    return ret.second;
  }
  
protected:

  std::unordered_set<smt::Term> walked_nodes_;
  static parent_map_t parent_;
  
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
public: // -- static -- model-to-model map
  using parent_map_t = ParentExtract::parent_map_t;

public:
  TermLearner(const smt::Term & trans_btor, to_next_t to_next_func, 
      //cex_term_map_t & cex_pred_map, 
      smt::SmtSolver & btor, const parent_map_t & parent_map) : 
    trans_(trans_btor), to_next_(to_next_func),
    // cex_pred_map_ref_(cex_pred_map),
    solver_(btor),
    parent_map_(parent_map) {}
    
  unsigned learn_terms_from_cex(IC3FormulaModel * pre, IC3FormulaModel * post, /*OUTPUT*/  PerVarsetInfo & varset_info );
  unsigned vars_extract_bit_level(IC3FormulaModel * post,  /*OUTPUT*/  PerVarsetInfo & varset_info) ;
  
protected:
  smt::Term trans_;
  to_next_t to_next_;
  //cex_term_map_t & cex_pred_map_ref_; // from unsat enum
  smt::SmtSolver & solver_;
  const parent_map_t & parent_map_; // ParentExtract
  
protected:
  unsigned same_val_replace_ast( /*INOUT*/  PerVarsetInfo & varset_info );
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

