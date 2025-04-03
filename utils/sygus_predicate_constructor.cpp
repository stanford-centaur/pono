/*********************                                                        */
/*! \file sygus_predicate_constructor.cpp
** \verbatim
** Top contributors (to current version):
**   Hongce Zhang
** This file is part of the pono project.
** Copyright (c) 2020 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Utility class to help manage/construct the predicates
**
**/

#include "utils/sygus_predicate_constructor.h"

#include <cassert>

namespace pono {
namespace syntax_analysis {

#define AND(x, y) (solver_->make_term(smt::And, (x), (y)))
#define OR(x, y) (solver_->make_term(smt::Or, (x), (y)))
#define NOT(x) (solver_->make_term(smt::Not, (x)))
#define EQ(x, y) (solver_->make_term(smt::Equal, (x), (y)))
#define LT(x, y) (solver_->make_term(smt::BVUlt, (x), (y)))
#define LE(x, y) (solver_->make_term(smt::BVUle, (x), (y)))
#define NEQ(x, y) (NOT(EQ((x), (y))))
#define TO_BV(x) (solver_->make_term(smt::Ite, (x), one, zero))

PredConstructor::PredConstructor(
    to_next_t to_next_func,
    smt::SmtSolver & solver,
    // const smt::Term & T_btor,
    // const smt::Term & Init_btor,
    // const smt::Term & Fprev_btor,
    IC3FormulaModel * cex,
    PerCexInfo & per_cex_info  // from pdr class, this will allow us to not use
                               // static data
    // bool init_per_cex_info,
    // VarTermManager & var_term_extractor
    // setup_cex_info will be in pdr class
    )
    : to_next_(to_next_func),
      solver_(solver),
      cex_(cex),
      per_cex_info_(per_cex_info),
      zero(solver_->make_term(0, solver_->make_sort(smt::SortKind::BV, 1))),
      one(solver->make_term(1, solver_->make_sort(smt::SortKind::BV, 1)))
{
  TermsDumping();
  terms_to_predicates();
}

smt::Term PredConstructor::smart_EQ(const smt::Term & l, const smt::Term & r)
{
  if (l->get_sort()->get_sort_kind() == smt::SortKind::BOOL
      && r->get_sort()->get_sort_kind() == smt::SortKind::BV
      && r->get_sort()->get_width() == 1) {
    return EQ(TO_BV(l), r);
  }
  if (r->get_sort()->get_sort_kind() == smt::SortKind::BOOL
      && l->get_sort()->get_sort_kind() == smt::SortKind::BV
      && l->get_sort()->get_width() == 1) {
    return EQ(l, TO_BV(r));
  }
  return EQ(l, r);
}

smt::Term PredConstructor::smart_NEQ(const smt::Term & l, const smt::Term & r)
{
  return NOT(smart_EQ(l, r));
}

smt::Term PredConstructor::smart_LT(const smt::Term & l, const smt::Term & r)
{
  if (l->get_sort()->get_sort_kind() == smt::SortKind::BOOL
      && r->get_sort()->get_sort_kind() == smt::SortKind::BV
      && r->get_sort()->get_width() == 1) {
    return LT(TO_BV(l), r);
  }
  if (r->get_sort()->get_sort_kind() == smt::SortKind::BOOL
      && l->get_sort()->get_sort_kind() == smt::SortKind::BV
      && l->get_sort()->get_width() == 1) {
    return LT(l, TO_BV(r));
  }
  return LT(l, r);
}

smt::Term PredConstructor::smart_LE(const smt::Term & l, const smt::Term & r)
{
  if (l->get_sort()->get_sort_kind() == smt::SortKind::BOOL
      && r->get_sort()->get_sort_kind() == smt::SortKind::BV
      && r->get_sort()->get_width() == 1) {
    return LE(TO_BV(l), r);
  }
  if (r->get_sort()->get_sort_kind() == smt::SortKind::BOOL
      && l->get_sort()->get_sort_kind() == smt::SortKind::BV
      && l->get_sort()->get_width() == 1) {
    return LE(l, TO_BV(r));
  }
  return LE(l, r);
}

uint64_t PredConstructor::term_summary() const
{
  return per_cex_info_.varset_info.terms_strings.size();
}

#define TERM_TABLE_DEBUG_LVL 0
void PredConstructor::TermsDumping() const
{
#if TERM_TABLE_DEBUG_LVL >= 1

  for (auto && w_term_cnst_pair : per_cex_info_.varset_info.terms) {
    auto width = w_term_cnst_pair.first;
    const PerWidthInfo & terms_consts = w_term_cnst_pair.second;

    auto nt_size = terms_consts.terms.size();
    auto nc_size = terms_consts.constants.size();

    std::cout << "[Width = " << width << "] " << "[#Term = " << nt_size
              << ", #Consts = " << nc_size << "]\n";

    std::cout << "  C : ";
#if TERM_TABLE_DEBUG_LVL < 2
    if (nc_size == 0)
      std::cout << " (none) " << std::endl;
    else if (nc_size == 1)
      std::cout << ((terms_consts.constants.at(0))->to_string()) << " (1) "
                << std::endl;
    else
      std::cout << (terms_consts.constants.front()->to_string()) << " .. "
                << (terms_consts.constants.back()->to_string()) << " ("
                << nc_size << ")" << std::endl;
#else
    for (auto && t : terms_consts.constants)
      std::cout << t->to_string() << ", ";
    std::cout << " (" << nc_size << ")\n";
#endif

    std::cout << "  T : ";
#if TERM_TABLE_DEBUG_LVL < 2
    if (nt_size == 0)
      std::cout << " (none) " << std::endl;
    else if (nt_size == 1)
      std::cout << ((terms_consts.terms.at(0))->to_string()) << " (1) "
                << std::endl;
    else
      std::cout << (terms_consts.terms.front()->to_string()) << " .. "
                << (terms_consts.terms.back()->to_string()) << " (" << nt_size
                << ")" << std::endl;
#else
    for (auto && t : terms_consts.terms) std::cout << t->to_string() << ", ";
    std::cout << " (" << nt_size << ")\n";
#endif
  }
#endif
}  // term dumping

#define ADD_PRED(pred_curr)                             \
  if (!((pred_curr)->is_value())) {                     \
    auto pred_next = to_next_(pred_curr);               \
    if (pred_str.insert(pred_curr->to_string()).second) \
      preds.push_back(pred_next);                       \
  }

// next_to_curr.emplace(pred_next, pred_curr);

void PredConstructor::terms_to_predicates()
{
  bool use_lt = per_cex_info_.varset_info.use_lt();
  bool use_lte = per_cex_info_.varset_info.use_lte();

  auto & preds = per_cex_info_.predicates_nxt;
  auto & pred_str = per_cex_info_.predicates_str;
  // auto & next_to_curr = per_cex_info_.pred_next_to_pred_curr;
  const auto & value_map = per_cex_info_.terms_val_under_cex;
  for (auto && w_term_cnst_pair : per_cex_info_.varset_info.terms) {
    // const - term
    auto width = w_term_cnst_pair.first;
    const PerWidthInfo & terms_consts = w_term_cnst_pair.second;
    const auto & terms = terms_consts.terms;
    const auto & constants = terms_consts.constants;

    unsigned nc = per_cex_info_.prev_per_width_term_num[width].const_num;
    unsigned nt = per_cex_info_.prev_per_width_term_num[width].term_num;
    auto nt_size = terms.size();
    auto nc_size = constants.size();

    for (unsigned cidx = 0; cidx < nc_size; ++cidx) {
      if (width <= 1 && cidx >= 1)
        break;  // you don't need a != 1 and a == 0, you only need to keep 1
                // value
      const auto & c = constants.at(cidx);
      const auto & cval = value_map.at(c);
      for (unsigned tidx =
               (cidx < nc
                    ? nt   // if c is an old one, we will start from new terms
                    : 0);  // else we can also use old terms
           tidx < nt_size;
           ++tidx) {
        const auto & t = terms.at(tidx);
        const auto & tval = value_map.at(t);

        auto pred_curr = (cval == tval) ? smart_EQ(c, t) : smart_NEQ(c, t);
        ADD_PRED(pred_curr)

        // use_lt
        if (use_lt && width > 1 && !(cval == tval)) {
          auto pred_curr = (cval < tval) ? smart_LT(c, t) : smart_LT(t, c);
          ADD_PRED(pred_curr)
        }

        if (use_lte && width > 1) {
          if (cval == tval || cval < tval) {
            auto pred = smart_LE(c, t);
            ADD_PRED(pred)
          }
          if (tval == cval || tval < cval) {
            auto pred = smart_LE(t, c);
            ADD_PRED(pred)
          }
        }  // use le
      }
    }  // end of c-t

    for (size_t idx1 = 0; idx1 < nt_size; ++idx1) {
      const auto & t1 = terms.at(idx1);
      const auto & tval1 = value_map.at(t1);
      for (size_t idx2 = (idx1 < nt ? nt : idx1 + 1);  // we must pick a new one
           idx2 < nt_size;
           ++idx2) {
        const auto & t2 = terms.at(idx2);
        const auto & tval2 = value_map.at(t2);
        auto pred_curr =
            (tval1 == tval2) ? smart_EQ(t1, t2) : smart_NEQ(t1, t2);
        ADD_PRED(pred_curr)

        // use_lt
        if (use_lt && width > 1 && !(tval1 == tval2)) {
          auto pred_curr =
              (tval1 < tval2) ? smart_LT(t1, t2) : smart_LT(t2, t1);
          ADD_PRED(pred_curr)
        }

        if (use_lte && width > 1) {
          if (tval1 == tval2 || tval1 < tval2) {
            auto pred = smart_LE(t1, t2);
            ADD_PRED(pred)
          }
          if (tval2 == tval1 || tval2 < tval1) {
            auto pred = smart_LE(t2, t1);
            ADD_PRED(pred)
          }
        }  // use le
      }
    }  // end of t-t

    // update record
    per_cex_info_.prev_per_width_term_num[width].const_num = nc_size;
    per_cex_info_.prev_per_width_term_num[width].term_num = nt_size;
  }  // end of w_term_cnst_pair

}  // terms_to_predicates

}  // namespace syntax_analysis
}  // namespace pono
