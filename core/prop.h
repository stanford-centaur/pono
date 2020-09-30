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

#include "smt-switch/smt.h"
#include "core/ts.h"

namespace pono {

class Property
{
 public:
  Property(TransitionSystem & ts, const smt::Term & p, std::string name="");

  /** Copies property to a new solver
   *  @param prop the property to copy
   *  @param tt the term translator to use
   *  @return a property using the solver in tt
   */
  Property(Property & prop, smt::TermTranslator & tt);

  ~Property();

  const smt::Term prop() const { return prop_; }

  TransitionSystem & transition_system() { return ts_; }

  std::string name() { return name_; };

 private:
  TransitionSystem ts_;

  const smt::Term prop_;

  std::string name_; ///< a name for the property. If no name is given, just uses the to_string

};  // class Property

}  // namespace pono
