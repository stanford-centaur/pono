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

#include "array_abstractor.h"

namespace pono {

class ArrayAbstractor : public Abstractor
{
 public:
  ArrayAbstractor(const TransitionSystem & ts, bool abstract_array_equality);

  smt::Term abstract(const smt::Term & t) const override;
  smt::Term concrete(const smt::Term & t) const override;

 protected:
  void do_abstraction() override;

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
