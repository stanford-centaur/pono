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

#include "core/prop.h"
#include "core/ts.h"
#include "core/unroller.h"
#include "modifiers/array_abstractor.h"
#include "refiners/axiom_enumerator.h"
#include "smt-switch/identity_walker.h"

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
  ARRAYEQ_READ_LAMBDA,
  LAMBDA_ALLDIFF
};

std::string to_string(AxiomClass ac);

// these are all the axioms that require instantiating an index
// crucial that this set is accurately maintained
// lambda axioms are not included because they're not parameterized
// by the index -- the index is known, lambda
// similarly, STORE_WRITE only uses the index in the store
const std::unordered_set<AxiomClass> index_axiom_classes(
    { CONSTARR, STORE_READ, ARRAYEQ_READ, LAMBDA_ALLDIFF });

// forward declaration for reference
class ArrayAxiomEnumerator;

// Walker for finding all the array terms and associated indices
// takes the *concrete* transition system and collects all array
// terms and indices and stores them in the appropriate
// data structures in the ArrayAxiomEnumerator
class ArrayFinder : public smt::IdentityWalker
{
  typedef smt::IdentityWalker super;

 public:
  ArrayFinder(ArrayAxiomEnumerator & aae);

 protected:
  smt::WalkerStepResult visit_term(smt::Term & term);

  ArrayAxiomEnumerator & aae_;
};

class ArrayAxiomEnumerator : public AxiomEnumerator
{
  friend ArrayFinder;

 public:
  ArrayAxiomEnumerator(ArrayAbstractor & aa,
                       Unroller & un,
                       const smt::Term & prop,
                       bool red_axioms);

  typedef AxiomEnumerator super;

  void initialize() override;

  bool enumerate_axioms(const smt::Term & abs_trace_formula,
                        size_t bound,
                        bool include_nonconsecutive = true) override;

  bool enumerate_axioms(const smt::Term & abs_trace_formula,
                        size_t bound,
                        bool include_nonconsecutive,
                        bool skip_lambda_axioms);

  /** Add a new index to the index set
   *  This can happen as we add auxiliary variables
   *  In particular, prophecy variables that are added here
   *  can help prove properties, even ones that require
   *  universally quantified invariants in some cases
   *
   *  @param the index to add to the index set
   */
  void add_index(const smt::Term & idx);

  /** Untimes an index, taking into account current and next
   *  e.g. the Unroller::untime(x@4 + y@5) would give x + y
   *  instead of x + y.next
   *  @param timed_index an unrolled index
   *  @return the untimed version
   */
  smt::Term untime_index(const smt::Term & timed_idx)
  {
    return untime_index_cache_.at(timed_idx);
  }

  smt::UnorderedTermSet & get_consecutive_axioms() override
  {
    return consecutive_axioms_;
  };

  AxiomVec & get_nonconsecutive_axioms() override
  {
    return nonconsecutive_axioms_;
  }

 protected:
  // helper functions

  /** populates all the data structures for generating axioms
   *  in a single traversal of the transition system
   */
  void collect_arrays_and_indices();

  /** creates lambda indices for each array index sort
   *  populates lambdas_
   */
  void create_lambda_indices();

  /** Clears all the data structures that are populated
   *  after a call to enumerate_axioms
   *  The expected use is that get_[non]consecutive_axioms
   *  only returns axioms from the latest abstract trace
   *  so it must be cleared before checking new axioms
   *  called in the beginning of enumerate_axioms
   */
  void clear_state();

  /** Check consecutive axioms from a certain class
   *  will populate consecutive_axioms_ with violated axioms
   *  @param ac the type of axiom to check
   *  @param only_curr if set to true then only checks axioms over current state
   * vars
   *  @param a limit on how many axioms to generate
   *         -1 means check all of them
   *  @return true iff any violated axioms were found
   */
  bool check_consecutive_axioms(AxiomClass ac,
                                bool only_curr,
                                int lemma_limit = -1);

  /** Check non-consecutive axioms from a certain class
   *  will populate nonconsecutive_axioms_ with violated axioms
   *  @param ac the type of axiom to check
   *  @param only_curr if set to true then only checks axioms over current state
   * vars
   *  @param i the time to instantiate indices at
   *  @param a limit on how many axioms to generate
   *         -1 means check all of them
   *  @return true iff any violated axioms were found
   */
  bool check_nonconsecutive_axioms(AxiomClass ac,
                                   bool only_curr,
                                   size_t i,
                                   int lemma_limit = -1);

  /** Check if a given axiom (over unrolled variables)
   *  is violated in the current model
   *  assumes the last call to the solver was satisfiable
   *  and there have been no pushes/pops since then
   *  @param ax the axiom to check
   *  @return true if the axiom is false in the current model
   */
  bool is_violated(const smt::Term & ax) const;

  // methods for instantiating groups of axioms
  // uses helper methods below for single axioms

  /** Instantiates axioms not in index_classes_
   *  i.e. they don't need a for loop over the index set
   *  @param ac the AxiomClass (assumed to not be in index_classes_)
   *  @return a set of axioms over transition system terms (not unrolled yet)
   */
  smt::UnorderedTermSet non_index_axioms(AxiomClass ac);

  /** Instantiates axioms in index_classes_
   *  i.e. will loop over indices
   *  @param ac the AxiomClass (assumed to not be in index_classes_)
   *  @param indices the set of indices to check (can be unrolled or not)
   *  @return a set of axioms over transition system terms (not - fully -
   * unrolled yet) Note: if checking non-consecutive axioms, the indices might
   * already be unrolled e.g. checking index i at a particular time
   */
  AxiomVec index_axioms(AxiomClass ac, smt::UnorderedTermSet & indices);

  // helper methods for instantiating single axioms

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

  /** Instantiates an all different axiom for lambdas
   *  from What's Decidable About Arrays
   *  the lambda index must be different from all indices of the same sort
   *  to handle the finite domain for bit-vectors, we actually use
   *  an integer sort for lambda, so that the all different constraint
   *  won't over-constrain (there's always a different integer)
   *  and convert bit-vectors to integers for comparison
   *  see lambda_guard
   *
   *  @param lambda the lambda index
   *  @param index another index from the index set (assumed to be of the same
   * concrete sort)
   *
   *  e.g. if lambda was created for the array sort (Array (_ BitVec 4) (_
   * BitVec 8)) then index should be the abstraction of an index of sort (_
   * BitVec 4)
   */
  smt::Term lambda_alldiff_axiom(const smt::Term & lambda,
                                 const smt::Term & index) const;

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

  /** Looks up the lambda for a particular index
   *  This is non-trivial because not all indices have a concrete
   *  version. For example, witnesses are only in the abstraction.
   *  Thus, they have no concrete sort to look up the lambda by
   *  but they're still associated with a particular array index sort
   *  this is tracked in witnesses_to_idxsort_
   *
   *  @param idx a non-lambda (abstract) index to find a corresponding lambda
   * for
   *  @return the lambda index corresponding to the same index sort
   */
  smt::Term get_lambda(smt::Term idx);

  /** Casts lambda to a bit-vector if necessary
   *  Currently bitvectors are the only finite-domain index supported
   *
   *  @param sort the concrete sort
   *  @param lam the lambda index
   *  @return the casted lambda so the sort matches the concrete sort
   */
  smt::Term cast_lambda(const smt::Sort & sort, const smt::Term & lam) const;

  // members
  // for abstracting/concretizing terms
  smt::Term conc_bad_;
  ArrayAbstractor & aa_;
  // for generating axioms
  Unroller & un_;

  bool reduce_axioms_unsatcore_;  ///< reduce generated axioms with an unsat
                                  ///< core if set to true

  size_t bound_;  ///< the bound of the current abstract trace
  smt::UnorderedTermMap
      constarrs_;  ///< maps (abstract) constarrs to their constant value
  smt::UnorderedTermSet stores_;  ///< vector of (abstract) stores
  // for index set, witness and lambda information
  // see What's Decidable About Arrays paper
  // the index set here does not contain lambdas
  // those need to be added separately for correctness
  smt::UnorderedTermSet index_set_;  ///< index set
  smt::UnorderedTermSet
      cur_index_set_;  ///< subset of index sets with terms containing only
                       ///< current state variables
  smt::UnorderedTermMap arrayeq_witnesses_;  ///< witnesses for array equalities
  std::unordered_map<smt::Term, smt::Sort>
      witnesses_to_idxsort_;  ///< maps witness index to corresponding concrete
                              ///< index sort
  std::unordered_map<smt::Sort, smt::Term>
      lambdas_;  ///< map from (concrete) array index sort to corresponding
                 ///< lambda

  // for axiom checking and storing
  smt::UnorderedTermSet
      axioms_to_check_;  ///< member variable used to store up axioms
  smt::UnorderedTermSet
      violated_axioms_;  ///< keeps track of violated axioms in given trace
  smt::UnorderedTermMap
      ts_axioms_;  ///< maps unrolled axioms to the transition system axioms
  std::unordered_map<smt::Term, AxiomInstantiation>
      to_axiom_inst_;  ///< maps a (violated) axiom term to its
                       ///< AxiomInstantiation object
  smt::UnorderedTermSet
      consecutive_axioms_;          ///< populated with consecutive axioms over
                                    ///< transition system variables
  AxiomVec nonconsecutive_axioms_;  ///< populated with nonconsecutive axiom
                                    ///< instantiations

  smt::UnorderedTermMap untime_index_cache_;

  // useful terms
  smt::Term false_;
};

}  // namespace pono
