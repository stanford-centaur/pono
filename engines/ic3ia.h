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
#include "modifiers/implicit_predicate_abstractor.h"
#include "smt-switch/term_translator.h"

namespace pono {

class IC3IA : public IC3
{
 public:
  // itp_se is the SolverEnum for the interpolator

  IC3IA(Property & p,
        smt::SolverEnum se,
        smt::SolverEnum itp_se = smt::SolverEnum::MSAT);
  IC3IA(Property & p,
        const smt::SmtSolver & s,
        smt::SolverEnum itp_se = smt::SolverEnum::MSAT);
  IC3IA(const PonoOptions & opt,
        Property & p,
        smt::SolverEnum se,
        smt::SolverEnum itp_se = smt::SolverEnum::MSAT);
  IC3IA(const PonoOptions & opt,
        Property & p,
        const smt::SmtSolver & s,
        smt::SolverEnum itp_se = smt::SolverEnum::MSAT);
  virtual ~IC3IA() {}

  typedef IC3 super;

 protected:
  // Note: important that conc_ts_ and abs_ts_ are before ia_
  //       because we will pass them to ia_ and they must be
  //       be initialized first

  TransitionSystem & conc_ts_;  ///< convenient reference to the concrete ts

  RelationalTransitionSystem abs_ts_;  ///< the abstract ts
                                       ///< after initialize, ts_ will point to
                                       ///< this because the methods from IC3
                                       ///< should operate on the abstraction

  ImplicitPredicateAbstractor ia_;

  smt::UnorderedTermSet predset_;  ///< set of current predicates
  // useful for checking if predicate has been added already
  // also available as a vector in ia_.predicates()

  smt::SmtSolver interpolator_;  ///< interpolator for refinement
  smt::TermTranslator
      to_interpolator_;  ///< transfer terms from solver_ to interpolator_
  smt::TermTranslator
      to_solver_;  ///< transfer terms from interpolator_ to solver_

  // TODO: since we're already using unroller_ over solver_ terms to reduce
  // predicates
  //       maybe better to just get rid of this and do all unrolling with
  //       solver_ terms should look into that
  TransitionSystem interp_ts_;  ///< ts_ over interpolator_ terms
  Unroller interp_unroller_;  ///< unroller for the interpolator for refinement

  // pure virtual method implementations

  IC3Formula get_ic3_formula(smt::TermVec * inputs = nullptr,
                             smt::TermVec * nexts = nullptr) const override;

  bool ic3_formula_check_valid(const IC3Formula & u) const override;

  // need to override this because IC3IA is not as restricted as
  // (bit-level) IC3
  void check_ts() const override;

  void initialize() override;

  void abstract() override;

  RefineResult refine() override;

  // specific to IC3IA

  /** Adds predicate to abstraction
   *  (calls ia_.add_predicate)
   *  and also incrementally updates the local transition relation
   *  and declares a new predicate state var (in pred_statevars_)
   *  @param pred the predicate over current state variables
   *  @return true iff the predicate was new (not seen before)
   */
  bool add_predicate(const smt::Term & pred);
};

}  // namespace pono
