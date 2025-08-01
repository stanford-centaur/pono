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

#pragma once

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>

#include "core/functional_unroller.h"
#include "core/prop.h"
#include "core/refineresult.h"
#include "core/ts.h"
#include "engines/ic3.h"
#include "options/options.h"
#include "smt-switch/smt.h"

namespace pono {

using EquivalenceClasses =
    std::unordered_map<smt::Sort,
                       std::unordered_map<smt::Term, smt::UnorderedTermSet>>;

using TypedTerms = std::unordered_map<smt::Sort, smt::UnorderedTermSet>;

class IC3SA : public IC3
{
 public:
  IC3SA(const SafetyProperty & p,
        const TransitionSystem & ts,
        const smt::SmtSolver & solver,
        PonoOptions opt = PonoOptions());

  virtual ~IC3SA() {}

  typedef IC3 super;

 protected:
  bool compute_witness() override;

  TransitionSystem conc_ts_;

  FunctionalUnroller f_unroller_;

  smt::UnorderedTermSet predset_;  ///< stores all predicates in abstraction

  // TODO remove this and generate it on the fly
  //      less performant from a time perspective but can get
  //      better generalization by using only terms in the actual query
  TypedTerms term_abstraction_;
  ///< stores all the current terms in the abstraction organized by sort

  smt::UnorderedTermSet projection_set_;  ///< variables always in projection

  smt::SmtSolver interpolator_;
  std::unique_ptr<smt::TermTranslator> to_interpolator_;
  std::unique_ptr<smt::TermTranslator> from_interpolator_;

  size_t longest_unroll_;

  // temporary data structure for conjunctive_assumptions
  smt::TermVec tmp_;

  // useful sort
  smt::Sort boolsort_;

  /** Overriding the method. This will return the concrete_ts_ because ts_ is an
   *  different (i.e., relational) view of the concrete_ts_.
   *  TODO this can probably be simplified later
   *       for CegarOpsUf, it has to be this way because otherwise
   *       initialization doesn't work correctly, and ts_ is replaced by
   *       a functional system when it should be relational
   */
  TransitionSystem & prover_interface_ts() override { return conc_ts_; };

  // virtual method implementations

  IC3Formula get_model_ic3formula() const override;

  bool ic3formula_check_valid(const IC3Formula & u) const override;

  void predecessor_generalization(size_t i,
                                  const smt::Term & c,
                                  IC3Formula & pred) override;

  void check_ts() const override;

  void abstract() override;

  RefineResult refine() override;

  void initialize() override;

  // IC3SA specific methods

  RefineResult ic3sa_refine_functional(smt::Term & learned_lemma);

  RefineResult ic3sa_refine_value(smt::Term & learned_lemma);

  /** Helper function to create labels for conjuncts of a constraint
   *  and assume the implication
   *  @param term the term to break into conjuncts and assume labels for
   *  @param used_lbls a set to keep track of used labels
   *  @param lbls the vector of labels to add to
   *  @param assumps the vector of assumptions to add to
   */
  void conjunctive_assumptions(const smt::Term & term,
                               smt::UnorderedTermSet & used_lbls,
                               smt::TermVec & lbls,
                               smt::TermVec & assumps);

  /** Get equivalence classes over all current terms in term_abstraction_
   *  from the current model
   *  @requires solver_ state is sat
   *  @param to_keep - set of symbols to include in the partition
   *  @param ec EquivalenceClasses to populate
   *  @return EquivalenceClass partition of the current term abstraction
   */
  void get_equivalence_classes_from_model(const smt::UnorderedTermSet & to_keep,
                                          EquivalenceClasses & ec) const;

  /** Generate the literals for the partition given by ec and add to cube
   *  @param ec the equivalence classes partition
   *  @param out_cube set of formulae to add to
   */
  void construct_partition(const EquivalenceClasses & ec,
                           smt::UnorderedTermSet & out_cube) const;

  /** Add all subterms from term to the term abstraction
   *  @param axiom the term to mine for subterms
   *  @return a set of new terms
   *  @modifies term_abstraction_ and predset_
   */
  smt::UnorderedTermSet add_to_term_abstraction(const smt::Term & term);

  /** Identifies a set of variables to project an abstract state onto
   *  based on justification (controlling arguments) and COI
   *  given the current model
   *  NOTE: will include next state symbols
   *  @requires solver_ context to be satisfiable
   *  @param term the term to traverse for relevant variables
   *  @param projection the projection set of variables to add to
   *  @ensures no calls to solver_
   *  @ensures assignments to the variables added to projection are
   *           sufficient for a satisfying assignment (e.g., they
   *           constitute a partial model)
   */
  void justify_coi(smt::Term term, smt::UnorderedTermSet & projection);

  bool is_controlled(smt::PrimOp po, const smt::Term & val) const;

  // helper function for justify_coi
  smt::Term get_controlling(smt::Term t) const;

  /** Check if a term is in the projection
   *  @param t the term to check
   *  @param to_keep the symbols for the projection
   *  @return true iff t only contains symbols from to_keep
   *          or the global projection_set_
   */
  inline bool in_projection(const smt::Term & t,
                            const smt::UnorderedTermSet & to_keep) const
  {
    // TODO can improve this with a specialized function which returns false
    // as soon as it finds a symbol outside of to_keep when traversing t

    smt::UnorderedTermSet free_vars;
    get_free_symbolic_consts(t, free_vars);
    bool in_projection = true;
    for (const auto & fv : free_vars) {
      if (to_keep.find(fv) == to_keep.end()
          && projection_set_.find(fv) == projection_set_.end()) {
        // this term contains a symbol not in the to_keep set
        in_projection = false;
        break;
      }
    }
    return in_projection;
  }

  /** Register a state variable mapping in to_solver_
   *  This is a bit ugly but it's needed because symbols aren't created in
   * to_solver_ so it needs the mapping from interpolator_ symbols to solver_
   * symbols
   *  TODO look into a cleaner solution
   *  @param i the unrolling for state variables
   *         makes sure not to repeat work
   */
  void register_symbol_mappings(size_t i);

  // debug methods
  void debug_print_equivalence_classes(EquivalenceClasses ec) const;
  void debug_print_syntax_abstraction() const;
  void debug_print_cex_trace(const std::vector<IC3Formula> & cex);
};

}  // namespace pono
