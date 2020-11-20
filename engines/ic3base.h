/*********************                                                  */
/*! \file ic3base.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan, Florian Lonsing
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract base class implementation of IC3 parameterized by
**        the unit used in frames, pre-image computation, and inductive
**        and predecessor generalization techniques.
**
**/
#pragma once

#include "engines/prover.h"

namespace pono {

class IC3Unit
{
 public:
  IC3Unit(const smt::TermVec & c) : children_(c), negated_(false) {}
  virtual ~IC3Unit() {}

  /** Get a single term representation
   *  Depends on the unit, e.g. a Disjunction unit would be
   *  be an OR of all the children.
   */
  virtual smt::Term get_term() const = 0;

  const smt::TermVec & get_children() const { return children_; };

  bool is_negated() const { return negated_; };

 protected:
  smt::TermVec children_;
  bool negated_;

  /** Check if this is a valid instance of this type of IC3Unit
   *  e.g. a Clause would make sure all the children are literals
   */
  virtual bool check_valid() const = 0;
};

typedef std::unique_ptr<IC3Unit> (*IC3UnitCreator)(const smt::TermVec & terms);

struct ProofGoal
{
  // based on open-source ic3ia ProofObligation
  std::unique_ptr<IC3Unit> target;
  size_t idx;
  std::unique_ptr<ProofGoal> next;

  ProofGoal(std::unique_ptr<IC3Unit> u, size_t i, std::unique_ptr<ProofGoal> n)
      : target(std::move(u)), idx(i), next(std::move(n))
  {
  }
};

class IC3Base : public Prover
{
 public:
  /** IC3Base constructors take the normal arguments for a prover
   *  + a function that can create an IC3Unit
   *  Depending on the derived class IC3 implementation, the exact
   *  type of IC3Unit will differ: e.g. Clause, Disjunction
   */
  IC3Base(Property & p, smt::SolverEnum se, IC3UnitCreator ic);
  IC3Base(Property & p, const smt::SmtSolver & s, IC3UnitCreator ic);
  IC3Base(const PonoOptions & opt,
          Property & p,
          smt::SolverEnum se,
          IC3UnitCreator ic);
  IC3Base(const PonoOptions & opt,
          Property & p,
          const smt::SmtSolver & s,
          IC3UnitCreator ic);

  typedef Prover super;

  void initialize() override;

 protected:
  ///< a function to create an IC3Unit
  IC3UnitCreator mk_unit;

  ///< the frames data structure.
  ///< a vector of the given Unit template
  ///< which changes depending on the implementation
  std::vector<std::vector<IC3Unit>> frames_;

  ///< stack of outstanding proof goals
  std::vector<ProofGoal> proof_goals_;

  // labels for activating assertions
  smt::Term init_label_;       ///< label to activate init
  smt::Term trans_label_;      ///< label to activate trans
  smt::TermVec frame_labels_;  ///< labels to activate frames
};

}  // namespace pono
