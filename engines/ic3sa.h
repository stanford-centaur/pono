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

#include "engines/ic3.h"
#include "utils/fcoi.h"

namespace pono {

using EquivalenceClasses =
    std::unordered_map<smt::Sort,
                       std::unordered_map<smt::Term, smt::UnorderedTermSet>>;

class IC3SA : public IC3
{
 public:
  IC3SA(const Property & p,
        const TransitionSystem & ts,
        const smt::SmtSolver & solver,
        PonoOptions opt = PonoOptions());

  virtual ~IC3SA() {}

  typedef IC3 super;

 protected:
  smt::UnorderedTermSet predset_;  ///< stores all predicates in abstraction

  std::unordered_map<smt::Sort, smt::UnorderedTermSet> term_abstraction_;
  ///< stores all the current terms in the abstraction organized by sort

  smt::UnorderedTermSet projection_set_;  ///< variables always in projection

  smt::UnorderedTermSet vars_in_bad_;  ///< variables occurring in bad

  // TODO: replace this with partial_model
  FunctionalConeOfInfluence fcoi_;

  std::vector<smt::UnorderedTermMap>
      inputvars_at_time_;  ///< unrolled variables
                           ///< only the unconstrained variables
                           ///< e.g. input variables and state vars
                           ///< with no next state update

  // virtual method implementations

  IC3Formula get_model_ic3formula(
      smt::TermVec * out_inputs = nullptr,
      smt::TermVec * out_nexts = nullptr) const override;

  bool ic3formula_check_valid(const IC3Formula & u) const override;

  IC3Formula generalize_predecessor(size_t i, const IC3Formula & c) override;

  void check_ts() const override;

  bool intersects_bad(IC3Formula & out) override;

  void initialize() override;

  RefineResult refine() override;

  // IC3SA specific methods

  /** Computes the symbolic post-image of
   *  pi /\ ci under the current model
   *  @requires current solver_ state is SAT
   *  @param i the current step
   *  @param subst the substitution map to use and update
   *  @param last_model_vals map to keep track of model assignments for inputs
   *         in current model (for convenience maps latest unrolling -- which
   *         didn't exist for this model to the "untimed" input variables)
   *  @return the symbolic post image given the current model at step i+1
   */
  smt::TermVec symbolic_post_image(size_t i,
                                   smt::UnorderedTermMap & subst,
                                   smt::UnorderedTermMap & last_model_vals);

  /** Create fresh symbolic constants for input variables
   *  and state variables with no next state at time i
   *  functions like the unroller, but only unrolls
   *  unconstrained variables
   *  @param i the time-step
   */
  void gen_inputvars_at_time(size_t i);

  /** Get equivalence classes over all current terms in term_abstraction_
   *  from the current model
   *  @requires solver_ state is sat
   *  @param to_keep - set of symbols to include in the partition
   *  @return EquivalenceClass partition of the current term abstraction
   */
  EquivalenceClasses get_equivalence_classes_from_model(
      const smt::UnorderedTermSet & to_keep) const;

  /** Generate the literals for the partition given by ec and add to cube
   *  @param ec the equivalence classes partition
   *  @param out_cube vector of formulae to add to
   */
  void construct_partition(const EquivalenceClasses & ec,
                           smt::TermVec & out_cube) const;

  /** Add all subterms from term to the term abstraction
   *  @param axiom the term to mine for subterms
   *  @return true iff new terms were added
   *  @modifies term_abstraction_ and predset_
   */
  bool add_to_term_abstraction(const smt::Term & term);

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

  // debug methods
  void debug_print_equivalence_classes(EquivalenceClasses ec) const;
  void debug_print_syntax_abstraction() const;
  void debug_print_cex_trace(const std::vector<IC3Formula> & cex);
};

}  // namespace pono
