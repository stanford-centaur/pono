/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#pragma once

#include "core/ts.h"
#include "smt-switch/smt.h"

namespace pono {

class Property
{
 public:
  Property(const smt::SmtSolver & s, const smt::Term & p, std::string name = "")
      : solver_(s), prop_vec_({ p }), name_(name)
  {
  }

  ~Property() {}

  Property(const smt::SmtSolver & s,
           const smt::TermVec & pv,
           std::string name = "")
      : solver_(s), prop_vec_(pv), name_(name)
  {
    if (pv.empty()) {
      throw PonoException("Property must have at least one term");
    }
  }

  const smt::Term & prop_term() const
  {
    if (prop_vec_.size() != 1) {
      throw PonoException("Property has multiple terms");
    }
    return prop_vec_.front();
  }

  const smt::TermVec & prop_vec() const { return prop_vec_; }

  const smt::SmtSolver & solver() const { return solver_; }

  std::string name() { return name_; };

 private:
  smt::SmtSolver solver_;

  smt::TermVec prop_vec_;

  std::string name_;  ///< a name for the property. If no name is given, just
                      ///< uses the to_string

};  // class Property

}  // namespace pono
