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
      : solver_(s), prop_(p), name_(name){};

  ~Property(){};

  const smt::Term & prop() const { return prop_; }

  const smt::SmtSolver & solver() const { return solver_; }

  std::string name() { return name_; };

 private:
  smt::SmtSolver solver_;

  smt::Term prop_;

  std::string name_;  ///< a name for the property. If no name is given, just
                      ///< uses the to_string

};  // class Property

}  // namespace pono
