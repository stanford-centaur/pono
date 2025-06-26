/*********************                                                  */
/*! \file axiom_enumerator.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract class for enumerating axioms over a transition system
**        it does not modify the transition system, instead only returning
**        violated axioms sufficient for ruling out abstract counterexamples.
**
**/
#pragma once

#include "core/ts.h"

namespace pono {

// Axiom instantiations
// a struct used to represent axioms that were instantiated with
// various terms
// what's considered an instantiation is up to the particular backend
// For Example,
//  the array theory has an axiom
//  forall a, b, i . a=b -> a[i] = b[i]
//  but the ArrayAxiomEnumerator only keeps track of the index (i)
//  instantiations and not the arrays a, b (because they don't need to be
//  tracked for later)
struct AxiomInstantiation
{
  AxiomInstantiation(const smt::Term & a, const smt::UnorderedTermSet & insts)
      : ax(a), instantiations(insts)
  {
  }

  smt::Term ax;  ///< the instantiated axiom
  smt::UnorderedTermSet
      instantiations;  ///< the instantiations used in the axioms
  // Note: the instantiations could be over unrolled variables, e.g. x@4 instead
  // of x
};

using AxiomVec = std::vector<AxiomInstantiation>;

class AxiomEnumerator
{
 public:
  // TODO: brainstorm the right interface
  //       would it be better for it to add axioms or not?
  AxiomEnumerator(const TransitionSystem & ts)
      : ts_(ts), solver_(ts.solver()), initialized_(false)
  {
  }

  virtual ~AxiomEnumerator(){};

  virtual void initialize() = 0;

  /** Check the axiom set over an abstract trace formula
   *  @param abs_trace_formula a formula representing the abstract trace
   *         it is assumed to be an (abstract) trace of ts_
   *  @param bound the bound the abstract trace was unrolled to
   *  @param include_nonconsecutive will look for non-consecutive axioms
   *         if set to true
   *         nonconsecutive means they cannot be "untimed" and added
   *         to the transition system directly
   *         e.g. a@10=b@10 -> select(a@10, i@4) = select(b@10, i@4)
   *         is non-consecutive because there's no way to add the axiom
   *         over current and next state variables because 10 - 4 > 1
   *  @return true iff the trace could be ruled out
   */
  virtual bool enumerate_axioms(const smt::Term & abs_trace_formula,
                                size_t bound,
                                bool include_nonconsecutive = true) = 0;

  /** Returns a sufficient set of violated consecutive axioms to
   *  rule out the abstract trace from the last call to
   *  enumerate_axioms
   *  @return a set of axiom instantiations that are consecutive
   *  meaning they only involve symbols from neighboring times
   *  and can be added directly to a transition system
   *  the free variables in the axiom instantiations are all
   *  state variables or inputs
   */
  virtual smt::UnorderedTermSet & get_consecutive_axioms() = 0;

  /** Returns a sufficient set of violated consecutive axioms to
   *  rule out the abstract trace from the last call to
   *  enumerate_axioms
   *  @return a vector of non-consecutive axiom instantiations
   *  These are structs that record extra information about the instantiation
   *  these axioms refer to timed symbols and cannot be added directly
   *  to a transition system.
   *  These must be handled with auxiliary variables or some other form
   *  of generalization such that the axiom becomes consecutive
   *  and can be added directly to the transition system.
   *  Examples:
   *
   *    a@4 = b@5 -> read(a@4, i@4) = read(b@5, i@4) is consecutive and can be
   *      added to the transition system as:
   *      a = b.next -> read(a, i) = read(b.next, i)
   *
   *    a@4 = b@5 -> read(a@4, i@7) = read(b@5, i@7) is non-consecutive
   *      because the mentioned times cannot be captured with just current
   *      and next.
   */
  virtual AxiomVec & get_nonconsecutive_axioms() = 0;

 protected:
  const TransitionSystem & ts_;
  const smt::SmtSolver & solver_;

  bool initialized_;
};

}  // namespace pono
