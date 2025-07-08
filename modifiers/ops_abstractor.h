/*********************                                                  */
/*! \file ops_abstractor.h
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan
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
#include "smt-switch/identity_walker.h"

namespace pono {

class OpsAbstractor : public Abstractor
{
 public:
  OpsAbstractor(const TransitionSystem & conc_ts, TransitionSystem & abs_ts);

  typedef Abstractor super;

  smt::Term abstract(smt::Term & t) override;
  smt::Term concrete(smt::Term & t) override;

  const smt::UnorderedTermMap & abstract_terms() const { return abs_terms_; }

  void set_ops_to_abstract(const smt::UnorderedOpSet & ops_to_abstract);

  void set_min_bitwidth(unsigned long w) { min_bw_ = w; };

  void do_abstraction();

 protected:
  class AbstractionWalker : smt::IdentityWalker
  {
   public:
    AbstractionWalker(OpsAbstractor & oa, smt::UnorderedTermMap * ext_cache);
    smt::Term visit(smt::Term & t) { return IdentityWalker::visit(t); }

   protected:
    smt::WalkerStepResult visit_term(smt::Term & t);
    OpsAbstractor & oa_;
  };
  friend class AbstractionWalker;

  class ConcretizationWalker : smt::IdentityWalker
  {
   public:
    ConcretizationWalker(OpsAbstractor & oa, smt::UnorderedTermMap * ext_cache);
    smt::Term visit(smt::Term & t) { return IdentityWalker::visit(t); }

   protected:
    smt::WalkerStepResult visit_term(smt::Term & t);
    OpsAbstractor & oa_;
  };
  friend class ConcretizationWalker;

  const smt::SmtSolver & solver_;

  AbstractionWalker abs_walker_;
  ConcretizationWalker conc_walker_;

  smt::UnorderedOpSet ops_to_abstract_;

  unsigned long min_bw_;

  std::unordered_map<std::string, smt::Term> abs_op_symbols_;
  std::unordered_map<smt::Term, smt::Op> abs_symbols_to_op_;

  smt::UnorderedTermMap abs_terms_;  // abs uf to concrete operator -- only
                                     // replace the top-level uf
};

}  // namespace pono
