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

#include "engines/prover.h"

namespace pono {

class KInduction : public Prover
{
 public:
  KInduction(const Property & p, smt::SolverEnum se);
  KInduction(const Property & p, const smt::SmtSolver & solver);
  KInduction(const PonoOptions & opt, const Property & p, smt::SolverEnum se);
  KInduction(const PonoOptions & opt,
             const Property & p,
             const smt::SmtSolver & solver);
  ~KInduction();

  typedef Prover super;

  void initialize() override;

  ProverResult check_until(int k) override;

 protected:
  bool base_step(int i);
  bool inductive_step(int i);

  smt::Term simple_path_constraint(int i, int j);
  bool check_simple_path_lazy(int i);

  smt::Term init0_;
  smt::Term false_;
  smt::Term simple_path_;

};  // class KInduction

}  // namespace pono
