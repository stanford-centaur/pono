/*********************                                                  */
/*! \file op_abstractor.h
** \verbatim
** Top contributors (to current version):
**   Hongce Zhang
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract complex operators with uninterpreted functions
**
**
**
**/

#pragma once

#include "core/unroller.h"
#include "engines/ic3base.h"
#include "modifiers/abstractor.h"
#include "smt-switch/identity_walker.h"
#include "smt-switch/utils.h"

namespace pono {

class OpAbstractor : public Abstractor
{
 public:
  OpAbstractor(const TransitionSystem & conc_ts, TransitionSystem & abs_ts)
      : Abstractor(conc_ts, abs_ts)
  {
  }

  // return true if it can refine
  // otherwise return false
  virtual bool refine_with_constraints(const smt::TermVec & cexs,
                                       const smt::Term & bad,
                                       smt::TermVec & out) = 0;

  virtual const smt::UnorderedTermSet & dummy_inputs() = 0;
  virtual bool has_abstracted() const = 0;
};  // class OpAbstractor

class OpInpAbstractor : public OpAbstractor
{
  struct OpAbstract
  {
    smt::Op op;
    smt::Term result;
    smt::TermVec args;

    smt::Term original;
    // statistics
    unsigned refine_count;

    OpAbstract() : refine_count(0) {}
  };

 public:
  typedef std::unordered_set<smt::PrimOp> OpSet;

  OpInpAbstractor(const TransitionSystem & conc_ts,
                  TransitionSystem & abs_ts,
                  const OpSet & op_to_abstract,
                  const smt::Term & prop,  // it is okay to use bad
                  int verbosity);

  // return true if it can refine
  // otherwise return false
  virtual bool refine_with_constraints(const smt::TermVec & cexs,
                                       const smt::Term & bad,
                                       smt::TermVec & out) override;

  virtual const smt::UnorderedTermSet & dummy_inputs() override
  {
    return dummy_inputs_;
  }

  virtual bool has_abstracted() const override
  {
    return !(op_abstracted.empty());
  }

 protected:
  // find mul/div and replace them with input/output

  smt::UnorderedTermSet dummy_inputs_;
  std::vector<OpAbstract> op_abstracted;

  void abstract_ts(const TransitionSystem & in_ts,
                   TransitionSystem & out_ts,
                   const OpSet & op_to_abstract,
                   const smt::Term & prop,  // it is okay to use bad
                   int verbosity);

  std::unique_ptr<Unroller> unroller_;

};  // class OpInpAbstractor

// ------------------------------------------------------
// extract uf in an expression

class UfExtractor : protected smt::IdentityWalker
{
 public:
  UfExtractor(const smt::SmtSolver & solver)
      : IdentityWalker(solver, true), out_(NULL)
  {
  }
  // need to clear the cache

  void extract_uf_in_term(const smt::Term & t, smt::UnorderedTermSet & uf_out);

 protected:
  smt::WalkerStepResult visit_term(smt::Term & term) override;

  smt::UnorderedTermSet * out_;
};  // class UfExtractor

class OpUfAbstractor : public OpAbstractor
{
  struct OpUfAbstract
  {
    smt::Op op;
    smt::Term result;
    smt::TermVec args;
    smt::Term uf;

    smt::Term original;
    // statistics
    unsigned refine_count;

    OpUfAbstract() : refine_count(0) {}
  };

 public:
  typedef std::unordered_set<smt::PrimOp> OpSet;

  OpUfAbstractor(const TransitionSystem & conc_ts,
                 TransitionSystem & abs_ts,
                 const OpSet & op_to_abstract);

  // return true if it can refine
  // otherwise return false
  virtual bool refine_with_constraints(const smt::TermVec & cexs,
                                       const smt::Term & bad,
                                       smt::TermVec & out) override;

  virtual const smt::UnorderedTermSet & dummy_inputs() override
  {
    return dummy_inputs_;
  }

  virtual bool has_abstracted() const override
  {
    return !(op_abstracted.empty());
  }

 protected:
  smt::UnorderedTermSet dummy_inputs_;  // empty set

  void abstract_ts(const TransitionSystem & in_ts,
                   TransitionSystem & out_ts,
                   const OpSet & op_to_abstract);

  std::unique_ptr<Unroller> unroller_;
  std::vector<OpUfAbstract> op_abstracted;
  std::unordered_map<std::string, smt::Term> uf_set_;
  UfExtractor uf_extractor_;
};  // class OpUfAbstractor

}  // namespace pono
