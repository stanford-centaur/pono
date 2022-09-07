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
  KInduction(const Property & p, const TransitionSystem & ts,
             const smt::SmtSolver & solver,
             PonoOptions opt = PonoOptions());

  virtual ~KInduction();

  typedef Prover super;

  void initialize() override;

  ProverResult check_until(int k) override;

 protected:
  smt::Term simple_path_constraint(int i, int j);
  bool check_simple_path_lazy(int i);
  bool check_simple_path_eager(int i);

  smt::Term init0_;
  smt::Term false_;
  smt::Term simple_path_;
  smt::Term neg_init_terms_;

  // Engine name used to print progress information used when running
  // k-induction and BMC + simple paths (a subclass of KInduction)
  std::string kind_engine_name_;
  // Wrapper function for logging feature. It takes the same
  // parameters as "logger.log(...)" except for string 'indent', which
  // is a string of space characters. The function calls
  // "logger.log(...)" and prepends the engine name to the output.
  template <typename... Args>
    void kind_log_msg(size_t level, const std::string & indent,
		      const std::string & format, const Args &... args);
};  // class KInduction

}  // namespace pono
