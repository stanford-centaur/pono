/*********************                                                  */
/*! \file ic3sa.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief IC3 with Syntax-Guided Abstraction based on
**
**        Model Checking of Verilog RTL Using IC3 with Syntax-Guided
*Abstraction.
**            -- Aman Goel, Karem Sakallah
**
**
**  within Pono, we are building on the bit-level IC3 instead of directly
**  on IC3Base, because a lot of the functionality is the same
**  In particular, we don't need to override inductive generalization
**
**/

#include "engines/ic3sa.h"

#include "assert.h"
#include "utils/term_walkers.h"

using namespace smt;
using namespace std;

namespace pono {

// main IC3SA implementation

IC3SA::IC3SA(Property & p, smt::SolverEnum se) : super(p, se) {}

IC3SA::IC3SA(Property & p, const smt::SmtSolver & s) : super(p, s) {}

IC3SA::IC3SA(const PonoOptions & opt, Property & p, smt::SolverEnum se)
    : super(opt, p, se)
{
}

IC3SA::IC3SA(const PonoOptions & opt, Property & p, const smt::SmtSolver & s)
    : super(opt, p, s)
{
}

IC3Formula IC3SA::get_model_ic3formula(TermVec * out_inputs,
                                       TermVec * out_nexts) const
{
  throw PonoException("IC3SA::get_model_ic3formula NYI");
}

bool IC3SA::ic3formula_check_valid(const IC3Formula & u) const
{
  throw PonoException("IC3SA::ic3formula_check_valid NYI");
}

IC3Formula IC3SA::generalize_predecessor(size_t i, const IC3Formula & c)
{
  throw PonoException("IC3SA::generalize_predecessor NYI");
}

void IC3SA::check_ts() const { throw PonoException("IC3SA::check_ts NYI"); }

RefineResult IC3SA::refine() { throw PonoException("IC3SA::refine NYI"); }

void IC3SA::initialize()
{
  super::initialize();

  // set up initial term abstraction by getting all subterms
  // TODO consider starting with only a subset -- e.g. variables
  SubTermCollector stc(solver_, false);
  stc.collect_subterms(ts_->init());
  stc.collect_subterms(ts_->trans());
  stc.collect_subterms(bad_);
  term_abstraction_ = stc.get_subterms();

  throw PonoException("IC3SA::initialize not completed");
}

// IC3SA specific methods

EquivalenceClasses IC3SA::get_equivalence_classes_from_model() const
{
  // assumes the solver state is sat
  EquivalenceClasses ec;
  for (auto elem : term_abstraction_) {
    const Sort & sort = elem.first;
    const UnorderedTermSet & terms = elem.second;

    // TODO figure out if a DisjointSet is a better data structure
    //      will need to keep track of all terms in each partition though

    ec[sort] = std::unordered_map<smt::Term, smt::UnorderedTermSet>();
    std::unordered_map<smt::Term, smt::UnorderedTermSet> & m = ec.at(sort);
    for (auto t : terms) {
      Term val = solver_->get_value(t);
      m[val].insert(t);
    }
  }
  return ec;
}

}  // namespace pono
