/*********************                                                        */
/*! \file sorting_network.h
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Sorting networks for symbolically sorting terms.
**
**
**/

#include "smt-switch/smt.h"

// based on implementation here:
// https://github.com/cristian-mattarei/CoSA/blob/master/cosa/utils/formula_mngm.py

namespace cosa {
// symbolically sorts a list of boolean terms
// true before false
class BoolSortingNetwork
{
 public:
  BoolSortingNetwork(smt::SmtSolver solver, smt::TermVec boolvec)
      : solver_(solver), boolvec_(boolvec)
  {
    do_sorting();
  }

  smt::TermVec sorted_boolvec() { return sorted_boolvec_; };

 protected:
  void do_sorting();
  smt::TermVec sorting_network_helper(smt::TermVec boolvec) const;
  smt::TermVec two_comparator(smt::Term bool0, smt::Term bool1) const;
  smt::TermVec merge(smt::TermVec boolvec0, smt::TermVec boolvec1) const;

  smt::SmtSolver solver_;
  smt::TermVec sorted_boolvec_;
  smt::TermVec boolvec_;
};
}  // namespace cosa
