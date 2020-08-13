/*********************                                                  */
/*! \file array_axiom_enumerator.cpp
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Class for enumerating array axioms over an array abstraction
**        produced by ArrayAbstractor (see array_abstractor.[h, cpp])
**
**/

#include "refiners/array_axiom_enumerator.h"

using namespace smt;
using namespace std;

namespace pono {

ArrayAxiomEnumerator::ArrayAxiomEnumerator(const TransitionSystem & ts,
                                           ArrayAbstractor & aa)
    : super(ts), aa_(aa)
{
}

bool ArrayAxiomEnumerator::enumerate_axioms(const Term & abs_trace_formula,
                                            size_t bound)
{
  // clear the previous axioms
  axioms_to_check_.clear();
  violated_axioms_.clear();
  ts_axioms_.clear();
  consecutive_axioms_.clear();
  nonconsecutive_axioms_.clear();

  solver_->push();
  throw PonoException("NYI");
  solver_->pop();
}

// protected methods

void ArrayAxiomEnumerator::collect_arrays_and_indices()
{
  throw PonoException("NYI");
}

void ArrayAxiomEnumerator::check_axioms(AxiomClass ac, int lemma_limit)
{
  throw PonoException("NYI");
}

bool ArrayAxiomEnumerator::is_violated(const Term & ax) const
{
  throw PonoException("NYI");
}

Term ArrayAxiomEnumerator::constarr_axiom(const Term & constarr,
                                          const Term & val,
                                          const Term & index) const
{
  throw PonoException("NYI");
}

Term ArrayAxiomEnumerator::constarr_lambda_axiom(const Term & constarr,
                                                 const Term & val) const
{
  throw PonoException("NYI");
}

Term ArrayAxiomEnumerator::store_write_axiom(const Term & store) const
{
  throw PonoException("NYI");
}

Term ArrayAxiomEnumerator::store_read_axiom(const Term & store,
                                            const Term & index) const
{
  throw PonoException("NYI");
}

Term ArrayAxiomEnumerator::store_read_lambda_axiom(const Term & store) const
{
  throw PonoException("NYI");
}

Term ArrayAxiomEnumerator::arrayeq_witness_axiom(const Term & arrayeq) const
{
  throw PonoException("NYI");
}

Term ArrayAxiomEnumerator::arrayeq_read_axiom(const Term & arrayeq,
                                              const Term & index) const
{
  throw PonoException("NYI");
}

Term ArrayAxiomEnumerator::arrayeq_read_lambda_axiom(const Term & arrayeq) const
{
  throw PonoException("NYI");
}

Term ArrayAxiomEnumerator::lambda_guard(const Sort & sort,
                                        const Term & lam) const
{
  throw PonoException("NYI");
}

}  // namespace pono
