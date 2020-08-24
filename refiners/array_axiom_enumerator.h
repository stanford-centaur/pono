/*********************                                                  */
/*! \file array_axiom_enumerator.h
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
#pragma once

#include "core/ts.h"
#include "modifiers/array_abstractor.h"
#include "refiners/axiom_enumerator.h"

namespace pono {

// enums representing each of the possible array axioms
enum AxiomClass
{
  CONSTARR = 0,
  CONSTARR_LAMBDA,
  STORE_WRITE,
  STORE_READ,
  STORE_READ_LAMBDA,
  ARRAYEQ_WITNESS,
  ARRAYEQ_READ,
  ARRAYEQ_READ_LAMBDA
};

class ArrayAxiomEnumerator : public AxiomEnumerator
{
 public:
  ArrayAxiomEnumerator(const TransitionSystem & ts, ArrayAbstractor & aa);

  typedef AxiomEnumerator super;

  bool enumerate_axioms(const smt::Term & abs_trace_formula,
                        size_t bound) override;

  smt::TermVec & get_consecutive_axioms() override
  {
    return consecutive_axioms_;
  };

  std::vector<NCAxiomInstantiation> & get_nonconsecutive_axioms() override
  {
    return nonconsecutive_axioms_;
  }

 protected:
  // helper functions

  /** populates all the data structures for generating axioms
   *  in a single traversal of the transition system
   */
  void collect_arrays_and_indices();

  /** Check axioms from a certain class
   *  @param ac the type of axiom to check
   *  @param a limit on how many axioms to generate
   *         -1 means check all of them
   */
  void check_axioms(AxiomClass ac, int lemma_limit = -1);

  /** Check if a given axiom (over unrolled variables)
   *  is violated in the current model
   *  assumes the last call to the solver was satisfiable
   *  and there have been no pushes/pops since then
   *  @param ax the axiom to check
   *  @return true if the axiom is false in the current model
   */
  bool is_violated(const smt::Term & ax) const;

  /** Instantiates the axiom:
   *
   *  forall i . select(constarr(val), i) = val
   *  at the given index over the abstracted constant array
   *
   *  @param constarr the abstract constant array
   *  @param val the element value of the concrete constant array
   *  @param index the index to instantiate the axiom at (can be unrolled or
   * not)
   *  @return the instantiated axiom
   */
  smt::Term constarr_axiom(const smt::Term & constarr,
                           const smt::Term & val,
                           const smt::Term & index) const;

  /** Instantiates the axiom:
   *
   *  forall i . select(constarr(val), i) = val
   *  at the lambda index
   *  it is very careful to guard the axiom appropriately if the
   *  (concrete) array sort has a finite domain index
   *  this is to avoid overconstraining issues where the entire
   *  domain is enumerated.
   *
   *  @param constarr the abstract constant array
   *  @param val the element value of the concrete constant array
   *  @return the instantiated axiom
   */
  smt::Term constarr_lambda_axiom(const smt::Term & constarr,
                                  const smt::Term & val) const;

  /** Creates the axiom:
   *
   *  select(store(a, j, e), j) = e
   *
   *  @param the abstract store
   *  @return the axiom
   */
  smt::Term store_write_axiom(const smt::Term & store) const;

  /** Instantiates the axiom:
   *
   *  forall i . i != j -> (select(store(a, j, e), i) = select(a, i))
   *  at the given index over the abstract arrays
   *
   *  @param store the abstract store term
   *  @param the index to instantiate it at (can be unrolled or not)
   *  @return the instantiated axiom
   */
  smt::Term store_read_axiom(const smt::Term & store,
                             const smt::Term & index) const;

  /** Instantiates the axiom:
   *
   *  forall i . i != j -> (select(store(a, j, e), i) = select(a, i))
   *  at the lambda index
   *  it is very careful to guard the axiom appropriately if the
   *  (concrete) array sort has a finite domain index
   *  this is to avoid overconstraining issues where the entire
   *  domain is enumerated.
   *
   *  @param store the abstract store term
   *  @return the instantiated axiom
   */
  smt::Term store_read_lambda_axiom(const smt::Term & store) const;

  /** Creates the axiom:
   *
   *  (a[witnesss] = b[witness]) -> a=b
   *    This is the only axiom that forces the
   *    arrays to be equal. Formally it's obtained from this lemma:
   *    (forall i . a[i] = b[i]) -> a = b
   *    !(forall i . a[i] = b[i]) | a = b
   *    (exists i . a[i] != b[i]) | a = b
   *    existential instantiation i -> witness
   *    a[witness] != b[witness] | a = b
   *    a[witness] = b[witness] -> a =b
   *
   *  @param arrayeq the abstract array equality (could be a UF depending on
   * options)
   *  @return the axiom
   */
  smt::Term arrayeq_witness_axiom(const smt::Term & arrayeq) const;

  /** Instantiates the axiom:
   *
   *  forall i . a = b -> a[i] = b[i]
   *  at the given index over the abstract arrays
   *
   *  @param arrayeq the abstract array equality
   *  @param the index to instantiate it at (can be unrolled or not)
   *  @return the instantiated axiom
   */
  smt::Term arrayeq_read_axiom(const smt::Term & arrayeq,
                               const smt::Term & index) const;

  /** Instantiates the axiom:
   *
   *  forall i . a = b -> a[i] = b[i]
   *  at the lambda index
   *  it is very careful to guard the axiom appropriately if the
   *  (concrete) array sort has a finite domain index
   *  this is to avoid overconstraining issues where the entire
   *  domain is enumerated.
   *
   *  @param arrayeq the abstract array equality
   *  @return the instantiated axiom
   */
  smt::Term arrayeq_read_lambda_axiom(const smt::Term & arrayeq) const;

  /** Creates the bounding guard for a lambda axiom
   *  for lambda's with an associated sort that has a
   *  finite domain. Currently should only be called with
   *  lambdas for bit-vector sorts
   *
   *  Example: if the index sort for this lambda is (_ BitVec 1)
   *  then there are only two possible values
   *  Thus, adding the constraint that lambda is different from all other
   * indices could make the queries trivially unsat Instead, we always use an
   * integer for the lambda And guard all lambda axioms with (0 <= lambda <=
   * upper_bound) -> axiom where in this case the upper bound is 1
   *
   *  @param sort the concrete array sort this lambda was created for
   *  @param lambda the lambda variable from the transition system
   *  @return a bounding guard for finite domains
   */
  smt::Term lambda_guard(const smt::Sort & sort, const smt::Term & lam) const;

  // members
  // for abstracting/concretizing terms
  ArrayAbstractor & aa_;
  // for generating axioms
  size_t bound_;  ///< the bound of the current abstract trace
  smt::UnorderedTermMap
      constarrs_;        ///< maps (abstract) constarrs to their constant value
  smt::TermVec stores_;  ///< vector of (abstract) stores
  smt::TermVec
      arrayeq_;  ///< vector of array equalities (could be UFs if abstracted)
  // for index set, witness and lambda information
  // see What's Decidable About Arrays paper
  // the index set here does not contain lambdas
  // those need to be added separately for correctness
  smt::UnorderedTermSet index_set_;  ///< index set
  smt::UnorderedTermSet
      cur_index_set_;  ///< subset of index sets with terms containing only
                       ///< current state variables
  smt::UnorderedTermMap arrayeq_witnesses_;  ///< witnesses for array equalities
  std::unordered_map<smt::Sort, smt::Term>
      lambdas_;  ///< map from (concrete) array sort to corresponding lambda

  // for axiom checking and storing
  smt::UnorderedTermSet
      axioms_to_check_;  ///< member variable used to store up axioms
  smt::UnorderedTermSet
      violated_axioms_;  ///< keeps track of violated axioms in given trace
  smt::UnorderedTermMap
      ts_axioms_;  ///< maps unrolled axioms to the transition system axioms
  smt::TermVec consecutive_axioms_;  ///< populated with consecutive axioms over
                                     ///< transition system variables
  std::vector<NCAxiomInstantiation>
      nonconsecutive_axioms_;  ///< populated with nonconsecutive axiom
                               ///< instantiations

  // useful terms
  smt::Term false_;
};

}  // namespace pono
