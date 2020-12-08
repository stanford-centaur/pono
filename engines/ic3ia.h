/*********************                                                  */
/*! \file ic3ia.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann, Ahmed Irfan
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief IC3 via Implicit Predicate Abstraction (IC3IA) implementation
**        based on
**
**        IC3 Modulo Theories via Implicit Predicate Abstraction
**            -- Alessandro Cimatti, Alberto Griggio,
**               Sergio Mover, Stefano Tonetta
**
**        and the open source implementation:
**
**        https://es-static.fbk.eu/people/griggio/ic3ia/index.html
**
**  within Pono, we are building on the bit-level IC3 instead of directly
**  on IC3Base, because a lot of the functionality is the same
**  In particular, we don't need to override either of the generalization
**  functions. Instead focusing on abstract/refine.
**
**/

#pragma once

#include "engines/ic3.h"

namespace pono {

class IC3IA : public IC3
{
 public:
  IC3IA(Property & p, smt::SolverEnum se);
  IC3IA(Property & p, const smt::SmtSolver & s);
  IC3IA(const PonoOptions & opt, Property & p, smt::SolverEnum se);
  IC3IA(const PonoOptions & opt, Property & p, const smt::SmtSolver & s);
  virtual ~IC3IA() {}

  typedef IC3 super;

 protected:
  // pure virtual method implementations

  IC3Formula get_ic3_formula() const override;

  bool ic3_formula_check_valid(const IC3Formula & u) const override;

  // need to override this because IC3IA is not as restricted as
  // (bit-level) IC3
  void check_ts() const override;

  void initialize() override;

  void abstract() override;

  RefineResult refine() override;
};

}  // namespace pono
