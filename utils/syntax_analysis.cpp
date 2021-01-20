/*********************                                                        */
/*! \file syntax_analysis.cpp
 ** \verbatim
 ** Top contributors (to current version):
 **   Hongce Zhang
 ** This file is part of the pono project.
 ** Copyright (c) 2020 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief Various functionalities for syntax analysis
 **
 ** 
 **/

#include "utils/syntax_analysis.h"
#include "utils/syntax_analysis_walker.h"
#include "utils/container_shortcut.h"

#include <cassert>

namespace pono {
namespace syntax_analysis {


// return delta(#terms)
// add options, add bool, add trans
unsigned VarTermManager::GetMoreTerms(IC3FormulaModel * pre, IC3FormulaModel * post, TermLearner & term_learner,
    const smt::Term & trans, bool failed_at_init, SyGuSTermMode term_mode) {
  // decide the policy
  assert(pre && post);
  std::string var_string = post->vars_to_canonical_string();
  assert(IN(var_string, terms_cache_)); // should already done this

  PerVarsetInfo & varset_info = terms_cache_.at(var_string);
  
  assert(varset_info.state.stage != PerVarsetInfo::state_t::EXTRACTBITS);
#if 0
  assert(varset_info.state.stage != PerVarsetInfo::state_t::VCLT);
  assert(varset_info.state.stage != PerVarsetInfo::state_t::VCLTE);
#endif
  assert(term_mode != SyGuSTermMode::VAR_C_EXT);
  // if we already extract bits, we should be forever good

  switch(varset_info.state.stage) {
    case PerVarsetInfo::state_t::VCLT:
    case PerVarsetInfo::state_t::VCLTE:
    case PerVarsetInfo::state_t::WPARTIAL:
      // increase to WALL
      {
        unsigned nterm_walked = insert_from_termsmap_w_width(
          varset_info.terms_buffer /*IN*/, varset_info /*OUT*/, varset_info.state.partial_width_done, (unsigned)(-1) );
        varset_info.state.stage = PerVarsetInfo::state_t::WALL;
        if (nterm_walked != 0)
          return nterm_walked;
      } // else will continue to WALL (same as from cex)
    case PerVarsetInfo::state_t::WALL:
      // increase to FROMCEX
      varset_info.state.stage = PerVarsetInfo::state_t::FROMCEX;
      // and do the same as FROMCEX, so continue
    case PerVarsetInfo::state_t::FROMCEX:
      // stay FROMCEX, just try if we can do anything more
      {
        auto nterms = 
        failed_at_init?
          term_learner.learn_terms_from_cex_same_frame(pre, post,
          trans, /*output*/ varset_info)
          :
          term_learner.learn_terms_from_cex(pre, post,
          trans, /*output*/ varset_info);
        if (nterms != 0)
          return nterms;
        // else
        varset_info.state.stage = PerVarsetInfo::state_t::EXTRACTBITS;
        // and will continue to the next stage
      }
    case PerVarsetInfo::state_t::EXTRACTBITS:
      return term_learner.vars_extract_bit_level(post, varset_info);
    default:
      assert(false);
  }
  assert(false);
} // GetMoreTerms

// we just won't compute canonical_string twice
const PerVarsetInfo & VarTermManager::SetupTermsForVarModelNormal(
  IC3FormulaModel * m, const std::string & canonical_string, 
  smt::SmtSolver & solver_,
  unsigned term_extract_depth, unsigned initial_term_width,
  unsigned initial_term_inc, unsigned accumulated_term_bound) {

  bool collect_constant = width_to_constants_.empty();
  std::unordered_set<smt::Term> varset;
  m->get_varset(varset);

  terms_cache_t::iterator pos; bool succ;
  std::tie(pos, succ) = terms_cache_.emplace(canonical_string, 
    PerVarsetInfo(PerVarsetInfo::state_t::EMPTY));

  // now TERM_EXTRACT_DEPTH
  TermExtractor extractor(varset, collect_constant, 
    term_extract_depth, 
    pos->second.terms_buffer,
    pos->second.all_terms);

  for (auto && t : terms_to_check_)
    extractor.WalkBFS(t);
  
  if (collect_constant)
    insert_from_constmap(extractor.GetConstants());
  // if collect constants

  auto & term_cache_item = pos->second;
  const auto & terms = pos->second.terms_buffer; // extractor.GetTerms();

  assert(!varset.empty());
  assert(!terms.empty());

  unsigned nterm_walked = 0;
  unsigned width_start = 0;
  unsigned width_end = initial_term_width;
  do{
    nterm_walked += insert_from_termsmap_w_width(terms /*IN*/, term_cache_item /*OUT*/, width_start, width_end );
    width_start += initial_term_inc;
    width_end += initial_term_inc;
  } while(nterm_walked <= accumulated_term_bound);

  pos->second.state.stage = PerVarsetInfo::state_t::WPARTIAL;
  pos->second.state.partial_width_done = width_start; // the start of next time

  term_const_w1_const(term_cache_item, solver_);

  return term_cache_item;
} // SetupTermsForVar


// will distinguish constants and terms because we don't need c==c or c=/=c
const PerVarsetInfo & VarTermManager::GetAllTermsForVarsInModel(
  IC3FormulaModel * m , smt::SmtSolver & s,
  SyGuSTermMode term_mode,
  unsigned term_extract_depth, unsigned initial_term_width,
  unsigned initial_term_inc, unsigned accumulated_term_bound) {

  std::string var_string = m->vars_to_canonical_string();
  auto pos = terms_cache_.find(var_string);
  if ( pos != terms_cache_.end() )  {
    return pos->second;
  }
  if (term_mode == SyGuSTermMode::FROM_DESIGN_LEARN_EXT)
    return SetupTermsForVarModelNormal(m, var_string, s, 
      term_extract_depth, initial_term_width, initial_term_inc, accumulated_term_bound);
  if (term_mode == SyGuSTermMode::VAR_C_EXT)
    return SetupTermsForVarModeExt(m, var_string, s);
  if (term_mode == SyGuSTermMode::SPLIT_FROM_DESIGN)
    return SetupTermsForVarModeSplit(m, var_string, s);
  if (term_mode == SyGuSTermMode::VAR_C_EQ_LT)
    return SetupTermsForVarModelVC(m, var_string, s); // just var and constant, you don't need a lot more

  assert(false);
} // GetAllTermsFor

// Later you may need to insert

// -------------------------------------------------------------------------


const PerVarsetInfo & VarTermManager::SetupTermsForVarModeSplit(
  IC3FormulaModel * m, const std::string & canonical_string, smt::SmtSolver & solver_)
{

  bool collect_constant = width_to_constants_.empty();
  std::unordered_set<smt::Term> varset;
  m->get_varset(varset);

  terms_cache_t::iterator pos; bool succ;
  std::tie(pos, succ) = terms_cache_.emplace(canonical_string, 
    PerVarsetInfo(PerVarsetInfo::state_t::WPARTIAL)); //WPARTIAL is better
    // EXTRACTBITS is simply see if we need more

  if (collect_constant) {
    // now TERM_EXTRACT_DEPTH
    ConstantExtractor extractor(width_to_constants_, constants_strings_);

    for (auto && t : terms_to_check_)
      extractor.WalkBFS(t);
    
    // -> width_to_constants_
    // todo make sure bit-1 has 1 and only 1 value
  }
  // if collect constants

  auto & term_cache_item = pos->second;
  // const auto & terms = term_cache_item.terms; // extractor.GetTerms();

  assert(!varset.empty());
  insert_vars_only(term_cache_item, varset);

  insert_split(term_cache_item, varset, solver_);
  unsigned nouse = 0;
  const_to_per_varset(term_cache_item, 0, (unsigned)(-1), nouse);

  // make sure there is at least a one for bv1
  term_const_w1_const(term_cache_item, solver_);

  return term_cache_item;
} // SetupTermsForVarModeSplit

const PerVarsetInfo & VarTermManager::SetupTermsForVarModelVC(
  IC3FormulaModel * m, const std::string & canonical_string,
  smt::SmtSolver & solver_)
{
  bool collect_constant = width_to_constants_.empty();
  std::unordered_set<smt::Term> varset;
  m->get_varset(varset);

  terms_cache_t::iterator pos; bool succ;
  std::tie(pos, succ) = terms_cache_.emplace(canonical_string, 
    PerVarsetInfo(PerVarsetInfo::state_t::VCLTE)); // see if we need more 
  // auto determine is needed!

  if (collect_constant) {
    // now TERM_EXTRACT_DEPTH
    ConstantExtractor extractor(width_to_constants_, constants_strings_);

    for (auto && t : terms_to_check_)
      extractor.WalkBFS(t);
    
    // -> width_to_constants_
    // todo make sure bit-1 has 1 and only 1 value
  }
  // if collect constants

  auto & term_cache_item = pos->second;
  // const auto & terms = term_cache_item.terms; // extractor.GetTerms();

  assert(!varset.empty());
  insert_vars_only(term_cache_item, varset);

  unsigned nouse = 0;
  const_to_per_varset(term_cache_item, 0, (unsigned)(-1), nouse);
  term_const_w1_const(term_cache_item, solver_);

  return term_cache_item;
} // SetupTermsForVarModelVC



const PerVarsetInfo & VarTermManager::SetupTermsForVarModeExt(
  IC3FormulaModel * m, const std::string & canonical_string, smt::SmtSolver & solver_) {
  
  bool collect_constant = width_to_constants_.empty();
  std::unordered_set<smt::Term> varset;
  m->get_varset(varset);

  terms_cache_t::iterator pos; bool succ;
  std::tie(pos, succ) = terms_cache_.emplace(canonical_string, 
    PerVarsetInfo(PerVarsetInfo::state_t::EXTRACTBITS));

  if (collect_constant) {
    // now TERM_EXTRACT_DEPTH
    ConstantExtractor extractor(width_to_constants_, constants_strings_);

    for (auto && t : terms_to_check_)
      extractor.WalkBFS(t);
    
    // -> width_to_constants_
    // todo make sure bit-1 has 1 and only 1 value
  }
  // if collect constants

  auto & term_cache_item = pos->second;
  // const auto & terms = term_cache_item.terms; // extractor.GetTerms();

  assert(!varset.empty());

  // insert vars and their extracts
  insert_vars_and_extracts(term_cache_item, varset, solver_);
  unsigned nouse = 0;
  const_to_per_varset(term_cache_item, 0, (unsigned)(-1), nouse);

  // make sure there is at least a one
  term_const_w1_const(term_cache_item, solver_);

  return term_cache_item;
} // SetupTermsForVar



// -------------------------------------------------------------------------

void VarTermManager::insert_split(
  PerVarsetInfo & term_cache_item /*OUT*/ , const smt::UnorderedTermSet & varset,
  smt::SmtSolver & solver_) 
{
  SliceExtractor extractor(term_cache_item.terms_buffer, varset);

  for (auto && t : terms_to_check_)
    extractor.WalkBFS(t);
  
  const auto & sv_slices = extractor.GetSvSlice();
  for (const auto & sv_slice_pair : sv_slices) {
    const auto & sv = sv_slice_pair.first;
    const auto & slices = sv_slice_pair.second;
    for (const auto & lr_pair : slices) {
      auto t = solver_->make_term(smt::Op(smt::PrimOp::Extract, lr_pair.first, lr_pair.second), sv);
      assert(lr_pair.first>=lr_pair.second);
      auto w = lr_pair.first-lr_pair.second+1;
      auto res = term_cache_item.terms_strings.insert(t->to_string());
      if (res.second)
        term_cache_item.terms[w].terms.push_back(t);
      // Parent relation? All_terms ?
    }
  }
} // insert_split


void VarTermManager::insert_vars_only(
  PerVarsetInfo & term_cache_item /*OUT*/ , const smt::UnorderedTermSet & varset)
{

  for (const auto & v : varset) {
    unsigned width;
    if (v->get_sort()->get_sort_kind() == smt::SortKind::BOOL )
      width = 1;
    else if (v->get_sort()->get_sort_kind() == smt::SortKind::BV )
      width = v->get_sort()->get_width();
    else
      continue;
    auto res = term_cache_item.terms_strings.insert(v->to_string());
    if(res.second)
      term_cache_item.terms[width].terms.push_back(v);
  }
} // insert_vars_only

void VarTermManager::insert_vars_and_extracts(
  PerVarsetInfo & term_cache_item /*OUT*/ , const smt::UnorderedTermSet & varset, smt::SmtSolver & solver_  ) 
{
  for (const auto & v : varset) {
    unsigned width;
    if (v->get_sort()->get_sort_kind() == smt::SortKind::BOOL )
      width = 1;
    else if (v->get_sort()->get_sort_kind() == smt::SortKind::BV )
      width = v->get_sort()->get_width();
    else
      continue;
    auto res = term_cache_item.terms_strings.insert(v->to_string());
    if(res.second) {
      term_cache_item.terms[width].terms.push_back(v);
    } // if insert successful    
    if (width > 1) {
      // make the extract
      for (unsigned idx = 0; idx < width; ++idx) {
        auto t = solver_->make_term(smt::Op(smt::PrimOp::Extract, idx, idx), v);
        auto res = term_cache_item.terms_strings.insert(t->to_string());
        if (res.second)
          term_cache_item.terms[1].terms.push_back(t);
      } // for each bit       
    } // if width > 1
  } // for each var
} // insert_vars_and_extracts

void VarTermManager::term_const_w1_const(PerVarsetInfo & term_cache_item /*OUT*/, smt::SmtSolver & solver_) {
  if (term_cache_item.terms.find(1) == term_cache_item.terms.end()  ||
        term_cache_item.terms.at(1).constants.empty()) {
    auto c0 = solver_->make_term(0);
    term_cache_item.terms_strings.insert(c0->to_string());
    term_cache_item.terms[1].constants.push_back(c0);
  }
} // term_const_w1_const

void VarTermManager::const_to_per_varset(PerVarsetInfo & term_cache_item /*OUT*/, 
  unsigned width_bound_low /*IN*/, unsigned width_bound_high /*IN*/, unsigned & nterm_walked)
{
  // insert constants to 
  for (auto && w_cs_pair : width_to_constants_) {
    auto width = w_cs_pair.first;
    if (! (width >= width_bound_low && width < width_bound_high))
      continue;

    const auto &cvec = w_cs_pair.second;
    for (auto && c : cvec) {
      auto ins_res = term_cache_item.terms_strings.insert(c->to_string());
      if(ins_res.second) {
        term_cache_item.terms[width].constants.push_back(c);
        ++ nterm_walked;
      } // end if really inserted
    }
  } // end of for each width_constant_pair
} //const_to_per_varset

unsigned VarTermManager::insert_from_termsmap_w_width(
  const std::map<unsigned, smt::TermVec> & terms /*IN*/, PerVarsetInfo & term_cache_item /*OUT*/ , 
  unsigned width_bound_low /*IN*/, unsigned width_bound_high /*IN*/) 
{
  unsigned nterm_walked = 0;
  for (auto && w_t_pair : terms) {
    auto width = w_t_pair.first;
    if (! (width >= width_bound_low && width < width_bound_high))
      continue;

    const auto & tvec = w_t_pair.second;
    for(auto && t : tvec) {
      auto tstr = t->to_string();
      auto ins_res = term_cache_item.terms_strings.insert(tstr);
      if (ins_res.second) { // if indeed inserted
        term_cache_item.terms[width].terms.push_back(t);
        ++ nterm_walked;
      }     
    } // end for each term
  } // end for terms
  const_to_per_varset(term_cache_item, width_bound_low, width_bound_high, nterm_walked);
  return nterm_walked;
} // insert_from_termsmap_w_width

void VarTermManager::insert_from_constmap(const PerVarsetInfo::width_term_map_t & w_c_map) {
  for (const auto & width_constvec_pair : w_c_map) {
    for (const auto & c : width_constvec_pair.second)  {
      auto cnstr_str = c->to_string();
      auto ins_res = constants_strings_.insert(cnstr_str);
      if (ins_res.second) // if insertion is successful
        width_to_constants_[width_constvec_pair.first].push_back(c);
    }
  }
} // insert_from_constmap

} // namespace syntax_analysis
} // namespace pono
