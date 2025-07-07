/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Ahmed Irfan, Makai Mann, Áron Ricardo Perez-Lopez
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief Classes representing LTL properties.
 **/

#pragma once

#include <string>

#include "smt-switch/smt.h"

namespace pono {

/**
 * \brief Abstract base class for all property types.
 *
 * Holds common metadata such as the associated solver and property name.
 */
class AbstractProperty
{
 public:
  /** @return the solver used to construct the property. */
  const smt::SmtSolver & solver() const;

  /** @return the name of this property. */
  const std::string name() const;

 protected:
  /** Construct a property. */
  AbstractProperty(const smt::SmtSolver & solver, const std::string name);

 private:
  const smt::SmtSolver solver_;

  const std::string name_;

};  // class AbstractProperty

/**
 * \brief Represents an LTL safety property.
 *
 * Corresponds to the negation of `bad` in Btor2.
 * A counterexample to a safety property is a finite path
 * from an initial state to a state where the property is violated
 * (i.e., a bad state).
 */
class SafetyProperty : public AbstractProperty
{
 public:
  /** Construct a safety property. */
  SafetyProperty(const smt::SmtSolver & solver,
                 const smt::Term & property,
                 std::string name = "");

  /** @return the SMT term of the property. */
  const smt::Term & prop() const;

 private:
  const smt::Term prop_;
};  // class SafetyProperty

using Property [[deprecated("Use SafetyProperty.")]] = SafetyProperty;

/**
 * \brief Represents an LTL liveness property.
 *
 * Corresponds to `justice` (and `fair`) in Btor2.
 * A counterexample to a liveness property is an
 * infinite lasso-shaped path starting from an initial state,
 * along which each term in the property is satisfied infinitely often.
 */
class LivenessProperty : public AbstractProperty
{
 public:
  /** Construct a liveness property. */
  LivenessProperty(const smt::SmtSolver & solver,
                   const smt::TermVec & conditions,
                   const std::string name = "");

  /** @return the terms (conditions) for the property. */
  const smt::TermVec & terms() const;

 private:
  const smt::TermVec prop_terms_;

};  // class LivenessProperty

}  // namespace pono
