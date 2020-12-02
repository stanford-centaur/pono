/*********************                                                  */
/*! \file ic3.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Bit-level IC3 implementation using the IC3Base abstract base class
**/

#pragma once

#include <unordered_map>
#include <unordered_set>

#include "engines/ic3base.h"

namespace pono {

class ClauseHandler : public IC3UnitHandler
{
 public:
  ClauseHandler(const smt::SmtSolver & s) : IC3UnitHandler(s) {}

  IC3Unit create(const smt::TermVec & c) const override;

  IC3Unit create_negated(const smt::TermVec & c) const override;

  IC3Unit negate(const IC3Unit & u) const override;

  bool check_valid(const IC3Unit & u) const override;
};

class IC3 : public IC3Base
{
 public:
  IC3(Property & p, smt::SolverEnum se);
  IC3(Property & p, const smt::SmtSolver & s);
  IC3(const PonoOptions & opt, Property & p, smt::SolverEnum se);
  IC3(const PonoOptions & opt, Property & p, const smt::SmtSolver & s);
  virtual ~IC3() {}

  typedef IC3Base super;

 protected:
  void initialize() override;

  // pure virtual method implementations

  std::vector<IC3Unit> inductive_generalization(size_t i,
                                                const IC3Unit & c) override;

  IC3Unit generalize_predecessor(size_t i, const IC3Unit & c) override;

  void check_ts() const override;

  IC3Unit get_unit() const override;
};

}  // namespace pono
