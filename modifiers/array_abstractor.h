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

#include "smt-switch/identity_walker.h"

#include "abstractor.h"

namespace pono {

class ArrayAbstractor;  // forward declaration for AbstractionWalker

class ArrayAbstractor : public Abstractor
{
  // Helper classes for walking formulas
  class AbstractionWalker : public smt::IdentityWalker
  {
  public:
    AbstractionWalker(ArrayAbstractor & aa, smt::UnorderedTermMap * ext_cache);
    smt::Term visit(smt::Term & t) { return IdentityWalker::visit(t); }

  protected:
    smt::WalkerStepResult visit_term(smt::Term & term);
    ArrayAbstractor & aa_;
  };

  class ConcretizationWalker : public smt::IdentityWalker
  {
  public:
    ConcretizationWalker(ArrayAbstractor & aa,
                         smt::UnorderedTermMap * ext_cache);
    smt::Term visit(smt::Term & t) { return IdentityWalker::visit(t); }

  protected:
    smt::WalkerStepResult visit_term(smt::Term & term);
    ArrayAbstractor & aa_;
  };

  friend class AbstractionWalker;
  friend class ConcretizationWalker;

 public:
  ArrayAbstractor(const TransitionSystem & conc_ts,
                  TransitionSystem & abs_ts,
                  bool abstract_array_equality = false);

  typedef Abstractor super;

  smt::Term abstract(smt::Term & t) override;
  smt::Term concrete(smt::Term & t) override;

  /** Returns the abstraction of a given sort
   *  if the sort has not been abstracted, the original
   *  sort is returned
   *  @param s the sort to look up the abstraction for
   *  @return the abstract sort
   */
  smt::Sort abstract(smt::Sort & s);
  /** Returns the concretization of a given abstract sort
   *  if the sort is not an abstraction, the original
   *  sort is returned
   *  @param s the sort to look up the concretization for
   *  @return the concrete sort
   */
  smt::Sort concrete(smt::Sort & s);

  /** Looks up a read UF for an abstract array sort
   *  @param sort the abstract array sort
   *  @return the corresponding read UF term
   */
  smt::Term get_read_uf(const smt::Sort & sort) const;
  /** Looks up a write UF for an abstract array sort
   *  @param sort the abstract array sort
   *  @return the corresponding write UF term
   */
  smt::Term get_write_uf(const smt::Sort & sort) const;
  /** Looks up a array equality UF for an abstract array sort
   *  @param sort the abstract array sort
   *  @return the corresponding array equality UF term
   *  NOTE: only valid if abstract_array_equality_ is true
   */
  smt::Term get_arrayeq_uf(const smt::Sort & sort) const;

  // getter
  // if true, then array equality is abstracted with a UF
  bool abstract_array_equality() const { return abstract_array_equality_; };

  void do_abstraction();

 protected:
  /** Populates abs_walker_'s cache and sort maps (abstract_sorts_,
   *  and concrete_sorts_) with array abstractions.
   *  Adds all the variables to the abs_ts_ --
   *    arrays are replaced with an abstraction.
   */
  void abstract_vars();

  /** Returns the abstract sort corresponding to a given
   *  array sort or creates one if necessary
   *  @param sort the concrete array sort
   *  @return the abstract array sort
   */
  smt::Sort abstract_array_sort(const smt::Sort & conc_sort);

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

  AbstractionWalker abs_walker_;
  ConcretizationWalker conc_walker_;

  ///< maps concrete sorts to abstract sorts
  std::unordered_map<smt::Sort, smt::Sort> abstract_sorts_;
  ///< maps abstract sorts to concrete sorts
  std::unordered_map<smt::Sort, smt::Sort> concrete_sorts_;

  ///< map from abstract array sort to read UF
  std::unordered_map<smt::Sort, smt::Term> read_ufs_;
  ///< map from abstract array sort to write UF
  std::unordered_map<smt::Sort, smt::Term> write_ufs_;
  ///< map from abstract array sort to equality UF
  std::unordered_map<smt::Sort, smt::Term> arrayeq_ufs_;

  // sets for the map values
  smt::UnorderedTermSet read_ufs_set_;
  smt::UnorderedTermSet write_ufs_set_;
  smt::UnorderedTermSet arrayeq_ufs_set_;
};

}  // namespace pono
