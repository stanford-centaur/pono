/*********************                                                  */
/*! \file static_coi.h
** \verbatim
** Top contributors (to current version):
**   Florian Lonsing, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Class for performing cone of influence reduction
**
**
**/

#pragma once

#include "core/ts.h"
#include "utils/fcoi.h"

namespace pono {
class StaticConeOfInfluence
{
 public:
  /** This class computes the cone of influence on construction
   *  @param ts the transition system to modify
   *  @param to_keep terms in the transition system that need to be kept
   *  The cone-of-influence will keep all the variables from terms in
   *    to_keep and any variables that influence those variables.
   */
  StaticConeOfInfluence(TransitionSystem & ts,
                        const smt::TermVec & to_keep,
                        int verbosity = 1);

 protected:

  TransitionSystem & ts_;
  int verbosity_;

  FunctionalConeOfInfluence coi_;  ///< class for computing symbols in
                                   ///< cone-of-influence of terms

  unsigned int orig_num_statevars_;
  unsigned int orig_num_inputvars_;
};
}  // namespace pono
