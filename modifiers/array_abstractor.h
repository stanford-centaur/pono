/*********************                                                  */
/*! \file array_abstractor.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract arrays using uninterpreted functions.
**
**
**/

#pragma once

#include "abstractor.h"

namespace pono {

class ArrayAbstractor : public Abstractor
{
 public:
  ArrayAbstractor(const TransitionSystem & ts,
                  bool abstract_array_equality = false);

  typedef Abstractor super;

  smt::Term abstract(const smt::Term & t) const override;
  smt::Term concrete(const smt::Term & t) const override;

 protected:
  void do_abstraction() override;

  /** Populates abstraction_cache_ and sort maps with array
   *  abstractions.
   */
  void abstract_array_vars();

  /** Returns the abstract sort corresponding to a given
   *  array sort or creates one if necessary
   *  @param sort the concrete array sort
   *  @return the abstract array sort
   */
  smt::Sort abstract_array_sort(const smt::Sort & conc_sort);

  /** Populates term caches from base class
   *  abstraction_cache_ and concretization_cache_
   *  @param conc_term the concrete term
   *  @param abs_term the abstract term
   *  asserts that values aren't overwritten
   */
  void update_term_cache(const smt::Term & conc_term,
                         const smt::Term & abs_term);

  /** Populates sort caches
   *  abstract_sorts_ and concrete_sorts_
   *  @param conc_sort the concrete sort
   *  @param abs_sort the abstract sort
   *  asserts that values aren't overwritten
   */
  void update_sort_cache(const smt::Sort & conc_sort,
                         const smt::Sort & abs_sort);

  bool abstract_array_equality_;
  const smt::SmtSolver & solver_;

  ///< maps concrete sorts to abstract sorts
  std::unordered_map<smt::Sort, smt::Sort> abstract_sorts_;
  ///< maps abstract sorts to concrete sorts
  std::unordered_map<smt::Sort, smt::Sort> concrete_sorts_;

  ///< map from abstract array sort to read UF
  std::unordered_map<smt::Sort, smt::Term> read_ufs_;
  ///< map from abstract array sort to write UF
  std::unordered_map<smt::Sort, smt::Term> write_ufs_;
  ///< map from abstract array sort to equality UF
  std::unordered_map<smt::Sort, smt::Term> eq_ufs_;
};

}  // namespace pono
