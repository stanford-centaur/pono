/*********************                                                        */
/*! \file syntax_analysis_walker.cpp
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

#include "utils/syntax_analysis_walker.h"
#include "utils/logger.h"
#include "utils/container_shortcut.h"
#include "utils/term_analysis.h"
 
#include "smt-switch/utils.h"

// #define DEBUG
#ifdef DEBUG
  #define D(...) logger.log( __VA_ARGS__ )
  #define INFO(...) D(0, __VA_ARGS__)
#else
  #define D(...) do {} while (0)
  #define INFO(...) logger.log(3, __VA_ARGS__)
#endif

namespace pono {
namespace syntax_analysis {
// I was worried that this will make stack overflow
void Walker::WalkDFS(const smt::Term & ast) {
  if (Skip(ast))
    return;
  PreChild(ast);
  for (const auto & c : *ast)
    WalkDFS(c);
  PostChild(ast);
}

void Walker::WalkBFS(const smt::Term & ast) {
  std::vector<std::pair<smt::Term,bool>> term_stack;
  
  term_stack.push_back(std::make_pair(ast,false));
  while(!term_stack.empty()) {
    auto & cur = term_stack.back();
    const auto & astnode = cur.first;
    if (Skip(astnode)) {
      term_stack.pop_back();
      continue;
    }

    if (cur.second) { // this is after its children
      PostChild(astnode);
      term_stack.pop_back();
      continue;
    } // this is after its child

    PreChild(astnode);
    cur.second = true;
    // for each of its child
    for (const auto & c : *astnode)
      term_stack.push_back(std::make_pair(c, false));
  }
} // end of Walker::WalkBFS

// -----------------------------------------------------


// ---------------------------------------------- //
//                                                //
//              OpExtractor                     //
//                                                //
// ---------------------------------------------- //

bool OpExtractor::Skip(const smt::Term & ast) {
  return IN(ast, walked_nodes_);
}

#define ARG1            \
      auto ptr = ast->begin();    \
      auto a1  = *(ptr++);      \
      assert (ptr == ast->end());

#define ARG2            \
      auto ptr = ast->begin();    \
      auto a1  = *(ptr++);      \
      auto a2  = *(ptr++);      \
      assert (ptr == ast->end());

void OpExtractor::PreChild(const smt::Term & ast) {
  walked_nodes_.insert(ast); // pre insert is okay
  auto sort_kind = ast->get_sort()->get_sort_kind() ;
  uint64_t width;
  
  if ( sort_kind == smt::SortKind::BOOL)
    width = 0;
  else if (sort_kind == smt::SortKind::BV)
    width = ast->get_sort()->get_width();
  else
    return; // we don't handle other cases


  if(ast->is_symbolic_const()) {
    // use ast->to_string()
    constructs.insert_symbol(width, ast->to_string());
  } else if (ast->is_value()) {
    // use ast->to_string()
    constructs.insert_const(width, ast->to_string());
  } else { // we will hope it is op
    smt::Op op;
    try {
      op = ast->get_op();
    } catch (NotImplementedException) {
      // op will be initialized correctly
      // so do nothing here
    }
  // op_unary, op_comp, op_concat, op_extract, op_rotate, op_extend
    // does not handle ite, repeat, convert, select, store ...
    switch (op.prim_op) {
      case smt::PrimOp::Not:
      case smt::PrimOp::BVNeg:
      case smt::PrimOp::BVNot: 
        constructs.insert_op_unary(width, op.prim_op); break;
      case smt::PrimOp::And:
      case smt::PrimOp::Or:
      case smt::PrimOp::Xor:
      case smt::PrimOp::Implies:
      // case smt::PrimOp::Iff:
      case smt::PrimOp::BVAnd:
      case smt::PrimOp::BVOr:
      case smt::PrimOp::BVXor:
      case smt::PrimOp::BVNor:
      case smt::PrimOp::BVXnor:
      case smt::PrimOp::BVNand:
      case smt::PrimOp::BVAdd:
      case smt::PrimOp::BVSub:
      case smt::PrimOp::BVMul:
      case smt::PrimOp::BVUdiv:
      case smt::PrimOp::BVSdiv:
      case smt::PrimOp::BVUrem:
      case smt::PrimOp::BVSrem:
      case smt::PrimOp::BVSmod:
      case smt::PrimOp::BVShl:
      case smt::PrimOp::BVAshr:
      case smt::PrimOp::BVLshr: 
        constructs.insert_op_binary(width, op.prim_op); break;
      case smt::PrimOp::BVComp: // equal 
        constructs.insert_op_comp(width, smt::PrimOp::Equal); break;
      case smt::PrimOp::Equal:
      case smt::PrimOp::Distinct:
      case smt::PrimOp::BVUlt:
      case smt::PrimOp::BVUle:
      case smt::PrimOp::BVUgt:
      case smt::PrimOp::BVUge:
      case smt::PrimOp::BVSlt:
      case smt::PrimOp::BVSle:
      case smt::PrimOp::BVSgt:
      case smt::PrimOp::BVSge: 
        constructs.insert_op_comp(width, op.prim_op); break;
      case smt::PrimOp::Concat:
        {
          ARG2
          constructs.insert_concat(width, concat_t(a1->get_sort()->get_width(),a2->get_sort()->get_width())); break;
        }
        break;
      case smt::PrimOp::Extract:
        assert (op.num_idx == 2);
        {
          ARG1
          constructs.insert_extract(width, extract_t(a1->get_sort()->get_width(),op.idx0,op.idx1)); 
        }
        break;
      
      case smt::PrimOp::Zero_Extend:
        assert (op.num_idx == 1);
        {
          ARG1
          constructs.insert_extend(width, extend_t(op.prim_op, op.idx0, a1->get_sort()->get_width()));
        }
        break;

      case smt::PrimOp::Sign_Extend:
        assert (op.num_idx == 1);
        {
          ARG1
          constructs.insert_extend(width, extend_t(op.prim_op, op.idx0, a1->get_sort()->get_width()));
        }
        break;


      case smt::PrimOp::Rotate_Left:
      case smt::PrimOp::Rotate_Right: {
        assert (op.num_idx == 1);
        constructs.insert_rotate(width, rotate_t(op.prim_op, op.idx0));
        break; }
      default: // do nothing
        break;
    } //  switch (op.prim_op)
  } // if it is op
} // OpExtractor::PreChild

void OpExtractor::PostChild(const smt::Term & ast) {

}

// ---------------------------------------------- //
//                                                //
//              TermExtractor                     //
//                                                //
// ---------------------------------------------- //

bool TermExtractor::Skip(const smt::Term & ast) {
  return IN(ast, walked_nodes_);
}

void TermExtractor::PreChild(const smt::Term & ast) {
  assert(!IN(ast, walked_nodes_));
}

void TermExtractor::PostChild(const smt::Term & ast) {
  // check if it is leaf

  walked_nodes_.emplace(ast, node_info_t() );

  unsigned width;
  auto sort_kind = ast->get_sort()->get_sort_kind() ;
  if ( sort_kind == smt::SortKind::BOOL)
    width = 1; // also make it bv?
  else if (sort_kind == smt::SortKind::BV)
    width = ast->get_sort()->get_width();
  else
    return ; // it is for example array, we don't handle

  if (ast->is_symbolic_const()) {
    // walked_nodes_[ast].level = 0; // default 0
    if (IN(ast,related_vars_)) {
      walked_nodes_[ast].in = true;
      walked_nodes_[ast].related = true;
      width_to_terms_[width].push_back(ast);
      all_terms_.insert(ast);
    }
  } else if ( ast->is_value() ) {

    walked_nodes_[ast].in = true;
    if (collect_constants_ ) {
      width_to_constants_[width].push_back(ast);
    }
  } else { // we will hope it is op
    unsigned max_level = 0;
    bool all_in = true;
    // bool some_in = false;
    D(0, "Walk : {} ", ast->to_string());

    for(auto && p : *ast) { // for each of its child node

      D(0, "  - Child : {} , lv: {} , in: {}", p->to_string(), walked_nodes_[p].level, walked_nodes_[p].in);
      max_level = std::max( walked_nodes_[p].level, max_level );
      all_in &= walked_nodes_[p].in;
    //  some_in |= walked_nodes_[p].related;
    }

    if (level_ > 0)
      ++ max_level; // if given level_ > 0, then we will count lv, otherwise, we don't care

    walked_nodes_[ast].in = all_in;
    // walked_nodes_[ast].related = some_in;
    walked_nodes_[ast].level = max_level;

    // D(0, "Result lv: {} , in: {} , related: {}", max_level, all_in, some_in);

    if (max_level <= level_ && all_in) {
      width_to_terms_[width].push_back(ast);
      all_terms_.insert(ast);
    }
    //if (max_level <= level_ && some_in) {
    //  related_terms_.insert(ast);
    //}
  } // end of op
} // PostChild


// ---------------------------------------------- //
//                                                //
//              TermScore Map                     //
//                                                //
// ---------------------------------------------- //


bool TermScore::Skip(const smt::Term & ast) {
  return IN(ast, scores_);
}

void TermScore::PreChild(const smt::Term & ast) {
  assert(!IN(ast, scores_));
 // walked_nodes_.insert(ast);
}

void TermScore::PostChild(const smt::Term & ast) {
 // for all its child, add parent pointer to the map
  unsigned width = 1;
  if (ast->get_sort()->get_sort_kind() == smt::SortKind::BOOL)
    width = 1;
  else if (ast->get_sort()->get_sort_kind() == smt::SortKind::BV)
    width = ast->get_sort()->get_width();

  if (ast->is_symbolic_const()) {
    scores_.emplace(ast,term_score_t(width)); // width*2
  } else if ( ast->is_value() ) { 
    scores_.emplace(ast,term_score_t(width*2)); // width
  } else { // we will hope it is op
    auto ret = scores_.emplace(ast,term_score_t(width));   // width  
    for(auto && c : *ast) { // for each of its child node
      ret.first->second.score += scores_.at(c).score;
      //  iterator->the score
    } // sum their scores and add one for itself
  } // end of op
} // PostChild

// ---------------------------------------------- //
//                                                //
//              Parent Extract                    //
//                                                //
// ---------------------------------------------- //

bool ParentExtract::Skip(const smt::Term & ast) {
  return IN(ast, walked_nodes_);
}

void ParentExtract::PreChild(const smt::Term & ast) {
  assert(!IN(ast, walked_nodes_));
 // walked_nodes_.insert(ast);
}

void ParentExtract::PostChild(const smt::Term & ast) {
 // for all its child, add parent pointer to the map
  walked_nodes_.insert(ast);
  if (ast->is_symbolic_const()) { }
  else if ( ast->is_value() ) { } 
  else { // we will hope it is op
    for(auto && p : *ast) { // for each of its child node
      parent_[p].insert(ast);
    } // set up its parent to have ast there
  } // end of op
} // PostChild

// ---------------------------------------------- //
//                                                //
//              Constant Extract                  //
//                                                //
// ---------------------------------------------- //

bool ConstantExtractor::Skip(const smt::Term & ast) {
  return (IN(ast, walked_nodes_));
}

void ConstantExtractor::PreChild(const smt::Term & ast) {
 // walked_nodes_.insert(ast);
}

void ConstantExtractor::PostChild(const smt::Term & ast) {
  walked_nodes_.insert(ast);

  if (! ast->is_value() )
    return;

  unsigned width;
  auto sort_kind = ast->get_sort()->get_sort_kind() ;
  if ( sort_kind == smt::SortKind::BOOL)
    width = 1; // also make it bv?
  else if (sort_kind == smt::SortKind::BV)
    width = ast->get_sort()->get_width();
  else
    return;
  auto ret = constants_strs_.insert(ast->to_string());
  if(ret.second)
    width_constant_map[width].push_back(ast);
} // PostChild

// ---------------------------------------------- //
//                                                //
//              SliceExtractor                    //
//                                                //
// ---------------------------------------------- //

bool SliceExtractor::Skip(const smt::Term & ast) {
  return IN(ast, walked_nodes_);
}

void SliceExtractor::PreChild(const smt::Term & ast) {

}

void SliceExtractor::PostChild(const smt::Term & ast) {
  walked_nodes_.insert(ast);

  if ( ast->is_value() )
    return;
  if ( ast->is_symbolic_const() )
    return;
  
  smt::Op op;
  try {
    op = ast->get_op();
  } catch (NotImplementedException) {
    // op will be initialized correctly
    // so do nothing here
    return;
  }

  if (op.prim_op == smt::PrimOp::Extract) {
    assert (op.num_idx == 2);

    // get its varset
    smt::UnorderedTermSet varset;
    smt::get_free_symbolic_consts(ast, varset);
    unsigned l = op.idx0, r = op.idx1;

    bool has_related_vars = false;
    for (const auto & sv : varset) {
      if (IN(sv, related_vars_)) {
        auto w = sv->get_sort()->get_width();
        if(l<w)
          sv2exts_[sv].insert(std::make_pair(l,r));
        has_related_vars = true;
      }
    } // for each var
    if (has_related_vars) {
      auto width = ast->get_sort()->get_width();
      //if (insert_res.second) { // if insertion succeed
      // this is just term buffer, does not matter
      ext_terms_[width].push_back(ast);
      //}
    }

  } // if it is extract (slice)
} // PostChild


// ---------------------------------------------- //
//                                                //
//                TermLearner                     //
//                                                //
// ---------------------------------------------- //



unsigned TermLearner::vars_extract_bit_level(IC3FormulaModel * post,  /*OUTPUT*/  PerVarsetInfo & varset_info) {
  smt::UnorderedTermSet vars;
  post->get_varset(vars);
  unsigned nterm = 0;
  for (const auto & v : vars) {
    if (!(v->get_sort()->get_sort_kind() == smt::SortKind::BV &&
          v->get_sort()->get_width() > 1 ))
      continue;
    auto width = v->get_sort()->get_width();
    for (unsigned idx = 0; idx < width; ++idx) {
      auto t = solver_->make_term(smt::Op(smt::PrimOp::Extract, idx, idx), v);
      parent_extractor_.RegisterNewParentRelation(v, t);
      bool is_new = varset_info.TermLearnerInsertTerm(t) ;
      //std::cout << "Extract " << v->to_string() << "[" << idx << "] is_new:" << is_new << std::endl;
      nterm += is_new ? 1 : 0;
    }
  }

  //GlobalTimer.RegisterEventCount("TermLearner.NewTermBl", nterm);
  return nterm;
} // vars_extract_bit_level



// return learned new terms
unsigned TermLearner::learn_terms_from_cex(
    IC3FormulaModel * pre_full_model, IC3FormulaModel * post,
    const smt::Term & trans,
    /*OUTPUT*/  PerVarsetInfo & varset_info ) {

  auto pre_prop = pre_full_model->to_expr();

  auto post_prop = solver_->make_term(smt::Not,(to_next_(post->to_expr())));
  unsigned delta_term_num = 0;
  D(0, "[TermLearner] Pre model : {}", pre_full_model->to_string() );
  D(0, "[TermLearner] Post model (will be negated) : {}", post->to_string() );
  solver_->push();
    solver_->assert_formula(pre_prop);
    solver_->assert_formula(trans);
    solver_->assert_formula(post_prop);
    auto res = solver_->check_sat();
    assert(res.is_sat());
    // okay now we need to find the right model on terms
    delta_term_num += concat_to_extract(varset_info);
    delta_term_num += same_val_replace_ast(varset_info); // only this use the model
    delta_term_num += extract_complement(varset_info);

  solver_->pop();
  D(0, "  [TermLearner] Learn new terms #{}", delta_term_num);
  return delta_term_num;  
} // learn_terms_from_cex


// return learned new terms
unsigned TermLearner::learn_terms_from_cex_same_frame(
    IC3FormulaModel * pre_full_model, IC3FormulaModel * post,
    const smt::Term & trans,
    /*OUTPUT*/  PerVarsetInfo & varset_info ) {

  auto pre_prop = pre_full_model->to_expr();

  auto not_post_prop = solver_->make_term(smt::Not,(post->to_expr()));
  unsigned delta_term_num = 0;
  D(0, "[TermLearner] Pre model : {}", pre_full_model->to_string() );
  D(0, "[TermLearner] Post model (will be negated) : {}", post->to_string() );
  solver_->push();
    solver_->assert_formula(pre_prop);
    solver_->assert_formula(not_post_prop);
    auto res = solver_->check_sat();
    assert(res.is_sat());
    // okay now we need to find the right model on terms
    delta_term_num += concat_to_extract(varset_info);
    delta_term_num += same_val_replace_ast_same_frame(varset_info); // only this use the model
    delta_term_num += extract_complement(varset_info);

  solver_->pop();
  D(0, "  [TermLearner] Learn new terms #{}", delta_term_num);
  return delta_term_num;  
} // learn_terms_from_cex

// helper function : get op
static smt::Op get_op(const smt::Term & ast) {
  smt::Op op;

  try {
    op = ast->get_op();
  } catch (NotImplementedException) {
    // op will be initialized correctly
    // so do nothing here
  }
  return op;
}

#undef ARG1
#undef ARG2
#define ARG1            \
      auto ptr = t->begin();    \
      auto a1  = *(ptr++);      \
      assert (ptr == t->end());

#define ARG2            \
      auto ptr = t->begin();    \
      auto a1  = *(ptr++);      \
      auto a2  = *(ptr++);      \
      assert (ptr == t->end());

unsigned TermLearner::concat_to_extract(/*INOUT*/  PerVarsetInfo & varset_info) {
  unsigned nterm = 0; 
  for(auto & width_info_pair : varset_info.terms){
    if (width_info_pair.first <= 1)
      continue; // you will not expect to have concat here
    const auto & terms = width_info_pair.second.terms;
    std::vector<std::pair<unsigned, unsigned>> extract_positions;
    for (const auto & t : terms) {
      // if it is concat, make extract ?
      if (!t->is_symbolic_const() && !t->is_value() && !t->is_param() &&
          t->begin() != t->end() && 
          (get_op(t).prim_op == smt::PrimOp::Concat)
          ) {
        
          ARG2
          unsigned sep = a2->get_sort()->get_width();
          unsigned msb = a1->get_sort()->get_width() + sep;
          // assert(sep>=1);
          extract_positions.push_back(std::make_pair(sep-1,0));
          extract_positions.push_back(std::make_pair(msb-1,sep));
      }
    } // end of for each term
    // okay, now, we will add the extract there
    for (const auto & t : terms) {
      for (const auto & pos : extract_positions) {
        auto new_term = solver_->make_term(smt::Op(smt::PrimOp::Extract, pos.first, pos.second), t);
        parent_extractor_.RegisterNewParentRelation(t, new_term);
        nterm += varset_info.TermLearnerInsertTerm(new_term) ? 1 : 0;
      } // for each postion
    } // for each term
  } // for each width
  return nterm;
} // concat_to_extract concat (v1, v2)  
// concat (extract(v, 6,0), extract(v,7,7) )
// concat (extract(v, 7,1), extract(v,0,0) )

unsigned TermLearner::extract_complement(/*INOUT*/  PerVarsetInfo & varset_info) {
  return 0; // not implemented
}

struct old_new_terms {
  smt::TermVec old;
  smt::TermVec noval;
};

#define REPLACE_DEBUG_LVL 0

#if REPLACE_DEBUG_LVL == 1
  #define RD1(...) D(__VA_ARGS__)
  #define RD2(...) do {} while(0)
#elif REPLACE_DEBUG_LVL == 2
  #define RD1(...) D(__VA_ARGS__)
  #define RD2(...) do {} while(0)
#else
  #define RD1(...) do {} while(0)  
  #define RD2(...) do {} while(0)  
#endif

unsigned TermLearner::same_val_replace_ast( /*INOUT*/  PerVarsetInfo & varset_info ) {
  unsigned n_new_terms = 0;
  for(auto & width_info_pair : varset_info.terms){
    auto width = width_info_pair.first;
    if (width <= 1)
      continue;
    std::unordered_map<eval_val, old_new_terms, eval_val_hash> eq_class; // equivalent classes

    for (const auto & t : width_info_pair.second.terms) {
      auto v_old = eval_val(solver_->get_value(t)->to_string());
      auto v_new = eval_val(solver_->get_value(to_next_(t))->to_string());
      eq_class[v_old].old.push_back(t);
      eq_class[v_new].noval.push_back(t);
    } // finish saparate into eq classes

    // also check for constants of the same width and insert
    for (auto const & c : width_info_pair.second.constants) {
      auto v = eval_val(c->to_string());
      eq_class[v].old.push_back(c);
      // eq_class[v].noval.push_back(c);
    }

    // for each pair, do a replacement and see if new terms can be added
    for( const auto & val_tvec_pair : eq_class) {
      const auto & tvec_old = val_tvec_pair.second.old;
      const auto & tvec_new = val_tvec_pair.second.noval;

#if REPLACE_DEBUG_LVL >= 1
      std::cout << "EQ class, val: " << val_tvec_pair.first.to_string() <<" w" << width
        << " #old:" << tvec_old.size() <<" |-> #new:" << tvec_new.size() << "\n  * ";
      for(const auto & t : tvec_old)
        std::cout <<t->to_string() << " , ";
      std::cout << " :|old -> new|: ";
      for(const auto & t : tvec_new)
        std::cout <<t->to_string() << " , ";
      std::cout << std::endl;
#endif
      // old -> new & new -> new
      for (const auto & oldterm : tvec_old) {
        for(const auto & newterm : tvec_new) {
          n_new_terms += replace_hierachically(oldterm, newterm, varset_info);
          if (n_new_terms > 1024)
            break;
        }
        if (n_new_terms > 1024)
          break;
      } // old -> new replacement

      // new -> new replacement
#if 0
      for( unsigned idx1 = 0; idx1 < tvec_new.size(); ++ idx1) {
        for (unsigned idx2 = idx1 + 1; idx2 < tvec_new.size(); ++ idx2) {
          const auto & t1 = tvec_new.at(idx1);
          const auto & t2 = tvec_new.at(idx2);
          
          n_new_terms += replace_hierachically(t1, t2, varset_info);
          n_new_terms += replace_hierachically(t2, t1, varset_info);
        }
      } // new -> new replacement


      for( unsigned idx1 = 0; idx1 < tvec_old.size(); ++ idx1) {
        for (unsigned idx2 = idx1 + 1; idx2 < tvec_old.size(); ++ idx2) {
          const auto & t1 = tvec_old.at(idx1);
          const auto & t2 = tvec_old.at(idx2);
          
          n_new_terms += replace_hierachically(t1, t2, varset_info);
          n_new_terms += replace_hierachically(t2, t1, varset_info);
        }
      } // new -> new replacement
#endif
    } // for each eq class
  } // on each width
  return n_new_terms;
} // same_val_replace_ast



unsigned TermLearner::same_val_replace_ast_same_frame( /*INOUT*/  PerVarsetInfo & varset_info ) {
  unsigned n_new_terms = 0;
  for(auto & width_info_pair : varset_info.terms){
    auto width = width_info_pair.first;
    if (width <= 1)
      continue;
    std::unordered_map<eval_val, old_new_terms, eval_val_hash> eq_class; // equivalent classes

    for (const auto & t : width_info_pair.second.terms) {
      auto v_old = eval_val(solver_->get_value(t)->to_string());
      eq_class[v_old].old.push_back(t);
    } // finish saparate into eq classes

    // also check for constants of the same width and insert
    for (auto const & c : width_info_pair.second.constants) {
      auto v = eval_val(c->to_string());
      eq_class[v].old.push_back(c);
    }

    // for each pair, do a replacement and see if new terms can be added
    for( const auto & val_tvec_pair : eq_class) {
      const auto & tvec_old = val_tvec_pair.second.old;

#if REPLACE_DEBUG_LVL >= 1
      std::cout << "EQ class, val: " << val_tvec_pair.first.to_string() <<" w" << width
        << " #old:" << tvec_old.size() <<" | replaced internally \n  * ";
      for(const auto & t : tvec_old)
        std::cout <<t->to_string() << " |end| ";
      std::cout << std::endl;
#endif

      // old -> old replacement
      for( unsigned idx1 = 0; idx1 < tvec_old.size(); ++ idx1) {
        for (unsigned idx2 = idx1 + 1; idx2 < tvec_old.size(); ++ idx2) {
          const auto & t1 = tvec_old.at(idx1);
          const auto & t2 = tvec_old.at(idx2);
          
          n_new_terms += replace_hierachically(t1, t2, varset_info);
          n_new_terms += replace_hierachically(t2, t1, varset_info);
        }
      } // old -> old replacement
      
    } // for each eq class
  } // on each width
  return n_new_terms;
} // same_val_replace_ast

unsigned TermLearner::replace_hierachically(
  const smt::Term & orig, const smt::Term & repl, PerVarsetInfo & varset_info ) {
  
  if (orig->to_string() == repl->to_string())
    return 0;

  smt::TermVec new_terms;
  RD1(3, "  [ReplaceInHierarchy] {} --> {} ", orig->to_string(), repl->to_string());
  unsigned orig_score = score_(orig), repl_score = score_(repl);
  RD1(3, "  [ReplaceInHierarchy] score {} --> {} ", orig_score, repl_score);
  if (repl_score >= orig_score * syntactic_score_factor + syntactic_score_delta) {
    RD1(3, "  [ReplaceInHierarchy] skipped. ");
    return 0;
  }

  unsigned ret = replace_hierachically_w_parent(orig, repl, varset_info, new_terms);
  for (const auto & nt : new_terms)
    RD1(3, "{}", nt->to_string() );
  RD1(3, "ret = {}", ret);
  assert(ret == new_terms.size());
  if (ret != 0) {
    for (const auto & t : new_terms) {
      parent_extractor_.WalkBFS(t);
    }
  }
  return ret;
}

// let's just make the stack larger
unsigned TermLearner::replace_hierachically_w_parent(
  const smt::Term & orig, const smt::Term & repl, PerVarsetInfo & varset_info,
  smt::TermVec & output_new_terms ) {
  
  RD2(3, "  [ReplaceParent] {} --> {} ", orig->to_string(), repl->to_string());
  const auto & parent_map_ = parent_extractor_.GetParentRelation();
  assert(! parent_map_.empty() );
  auto parent_termvec_pos = parent_map_.find(orig);
  if (parent_termvec_pos == parent_map_.end())
    return 0; // if has no parent term, no replacement
  const auto & parentvec = parent_termvec_pos->second;
  if (parentvec.empty())
    return 0; // if has no parent term, no replacement
  unsigned nterm = 0;

  for(const auto & p : parentvec ) {
    if (varset_info.TermLearnerIsOut(p)) {
      RD2(3, "  [ReplaceParent]    not in parent: {} , out", p->to_string() );
      continue;
    }

    smt::TermVec new_terms;
    {
      smt::TermVec old_children;
      std::vector<unsigned> child_pos;
      unsigned idx = 0;
      for(const auto & c : *p) {
        old_children.push_back(c);
        if (c == orig)
          child_pos.push_back(idx);
        ++ idx;
      } // find child pos
      for (auto idx : child_pos) {
        old_children[idx] = repl;
        auto new_parent = (solver_->make_term(p->get_op(), old_children));

       bool is_new_term = varset_info.TermLearnerInsertTerm(new_parent);
        if (is_new_term) {
          RD2(3, "  [ReplaceParent]    in parent (new): {} ==> {}", p->to_string(), new_parent->to_string() );
          output_new_terms.push_back(new_parent);
          new_terms.push_back(new_parent);
        } else {
          RD2(3, "  [ReplaceParent]    in parent (exists): {} ==> {}", p->to_string(), new_parent->to_string() );
        }
        old_children[idx] = orig;
      }
    } // replace_child_in_parent
    // avoid replace 1-0
    nterm += new_terms.size();
    if (nterm > 0xff)
      continue;
    if ( p->get_sort()->get_sort_kind() == smt::SortKind::BOOL ||
        (p->get_sort()->get_sort_kind() == smt::SortKind::BV &&
         p->get_sort()->get_width() == 1))
      continue;
    for (const auto & nt : new_terms) {
      RD2(1,"  [TermLearner Replace] {} ==> {}", p->to_string(), nt->to_string());
      nterm += replace_hierachically_w_parent(p, nt, varset_info,output_new_terms );
        //ParentExtract::RegisterNewParentRelation(c, out.back());
      if (nterm > 0xff)
        break;
    }
  }
  return nterm;
} // replace_hierachically

} // namespace syntax_analysis
} // namespace pono

