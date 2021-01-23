/*********************                                                        */
/*! \file ceg_prophecy_arrays.cpp
** \verbatim
** Top contributors (to current version):
**   Ahmed Irfan, Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief An implementation of Counter-Example Guided Prophecy for array
**        model checking. It is parameterized by an underlying model checking
**        procedure which need not handle arrays (only UF). However, a common
**        instantiation is with an IC3-style procedure, in which case we
**        often refer to this algorithm as "prophic3".
**
**/

#include "engines/ceg_prophecy_arrays.h"

#include <map>

#include "assert.h"
#include "engines/bmc.h"
#include "engines/bmc_simplepath.h"
#include "engines/ic3ia.h"
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"
#include "engines/mbic3.h"
#include "smt/available_solvers.h"

#include "utils/logger.h"
#include "utils/make_provers.h"
#include "utils/term_analysis.h"

#ifdef WITH_MSAT_IC3IA
#include "engines/msat_ic3ia.h"
#endif

using namespace smt;
using namespace std;

namespace pono {

template <class Prover_T>
CegProphecyArrays<Prover_T>::CegProphecyArrays(const Property & p,
                                               const TransitionSystem & ts,
                                               const SmtSolver & solver,
                                               PonoOptions opt)
    : super(p, RelationalTransitionSystem(solver), solver, opt),
      conc_ts_(ts),
      abs_ts_(super::prover_interface_ts()),
      abs_unroller_(abs_ts_, "{@}"),
      aa_(conc_ts_, abs_ts_, true),
      aae_(aa_,
           abs_unroller_,
           ts.solver() == super::solver_
               ? p.prop()
               : super::to_prover_solver_.transfer_term(p.prop(), BOOL),
           super::options_.cegp_axiom_red_),
      pm_(abs_ts_),
      reached_k_(-1),
      num_added_axioms_(0)
{
  // point orig_ts_ to the correct one
  super::orig_ts_ = ts;
}

template <class Prover_T>
ProverResult CegProphecyArrays<Prover_T>::prove()
{
  initialize();

  ProverResult res = ProverResult::FALSE;
  while (res == ProverResult::FALSE) {
    // Refine the system
    // heuristic -- stop refining when no new axioms are needed.
    do {
      if (!CegProphecyArrays::cegar_refine()) {
        // real counterexample
        return ProverResult::FALSE;
      }
      reached_k_++;
    } while (num_added_axioms_);

    if (super::engine_ != IC3IA_ENGINE) {
      Property latest_prop(super::solver_,
                           super::solver_->make_term(Not, super::bad_));
      SmtSolver s = create_solver_for(super::solver_->get_solver_enum(),
                                      super::engine_, false);
      shared_ptr<Prover> prover = make_prover(super::engine_, latest_prop,
                                              abs_ts_, s, super::options_);
      res = prover->prove();

      if (res == ProverResult::FALSE) {
        // use witness length
        // reached_k_ is the last k without a counterexample trace
        reached_k_ = prover->witness_length() - 1;
      }

    } else {
      res = super::prove();
      if (res == ProverResult::FALSE) {
        // use witness length
        reached_k_ = super::reached_k_;
      }
    }
  }

  if (res == ProverResult::TRUE && super::invar_) {
    // update the invariant
    super::invar_ = aa_.concrete(super::invar_);
  }

  return res;
}

template <class Prover_T>
ProverResult CegProphecyArrays<Prover_T>::check_until(int k)
{
  initialize();

  ProverResult res = ProverResult::FALSE;
  while (res == ProverResult::FALSE && reached_k_ <= k) {
    // Refine the system
    // heuristic -- stop refining when no new axioms are needed.
    do {
      if (!CegProphecyArrays::cegar_refine()) {
        return ProverResult::FALSE;
      }
      reached_k_++;
    } while (num_added_axioms_ && reached_k_ <= k);

    if (super::engine_ != IC3IA_ENGINE) {
      Property latest_prop(super::solver_,
                           super::solver_->make_term(Not, super::bad_));
      SmtSolver s = create_solver_for(super::solver_->get_solver_enum(),
                                      super::engine_, false);
      shared_ptr<Prover> prover = make_prover(super::engine_, latest_prop,
                                              abs_ts_, s, super::options_);
      res = prover->check_until(k);

      if (res == ProverResult::FALSE) {
        // use witness length
        // reached_k_ is the last k without a counterexample trace
        reached_k_ = prover->witness_length() - 1;
      } else if (res == ProverResult::TRUE) {
        try {
          // set the invariant
          super::invar_ = prover->invar();
        }
        catch (std::exception & e) {
          logger.log(3, "Failed to set invariant because {}", e.what());
          continue;
        }
      }
    } else {
      res = super::check_until(k);
      if (res == ProverResult::FALSE) {
        // use witness length
        reached_k_ = super::reached_k_;
      }
    }
  }

  if (res == ProverResult::FALSE) {
    // can't count on false result over abstraction when only checking up until
    // a bound
    return ProverResult::UNKNOWN;
  }

  if (res == ProverResult::TRUE && super::invar_) {
    // update the invariant
    super::invar_ = aa_.concrete(super::invar_);
  }

  return res;
}

template <class Prover_T>
void CegProphecyArrays<Prover_T>::initialize()
{
  if (super::initialized_) {
    return;
  }

  // specify which cegar_abstract in case
  // we're inheriting from another cegar algorithm
  CegProphecyArrays::cegar_abstract();
  // call super initializer after abstraction
  super::initialize();

  bool contains_arrays = false;
  for (const auto &sv : conc_ts_.statevars()) {
    if (sv->get_sort()->get_sort_kind() == ARRAY) {
      contains_arrays = true;
      break;
    }
  }

  for (const auto &iv : conc_ts_.inputvars()) {
    if (iv->get_sort()->get_sort_kind() == ARRAY) {
      contains_arrays = true;
      break;
    }
  }

  if (!contains_arrays) {
    throw PonoException("Ran CegProphecyArrays on system without arrays.");
  }
}

template <class Prover_T>
void CegProphecyArrays<Prover_T>::cegar_abstract()
{
  aa_.do_abstraction();
  aae_.initialize();
  // the ArrayAbstractor already abstracted the transition system on
  // construction -- only need to abstract bad
  assert(super::bad_);
  super::bad_ = aa_.abstract(super::bad_);

  // the abstract system should have all the same state and input variables
  // but abstracted
  // plus it will have some new variables for the witnesses and lambdas
  assert(abs_ts_.statevars().size() >= conc_ts_.statevars().size());
  assert(abs_ts_.inputvars().size() >= conc_ts_.inputvars().size());
}

template <class Prover_T>
bool CegProphecyArrays<Prover_T>::cegar_refine()
{
  num_added_axioms_ = 0;
  // TODO use ArrayAxiomEnumerator and modifiers to refine the system
  // create BMC formula
  Term abs_bmc_formula = get_bmc_formula(reached_k_ + 1);

  // check array axioms over the abstract system
  if (!aae_.enumerate_axioms(abs_bmc_formula, reached_k_ + 1)) {
    // concrete CEX
    return false;
  }

  UnorderedTermSet consecutive_axioms = aae_.get_consecutive_axioms();
  AxiomVec nonconsecutive_axioms = aae_.get_nonconsecutive_axioms();

  bool found_nonconsecutive_axioms = nonconsecutive_axioms.size();

  if (found_nonconsecutive_axioms) {
    // TODO: consider adding these axioms directly over the prophecy
    //       variables at the correct time
    //       for now, easier to just search for consecutive axioms

    if (super::options_.cegp_axiom_red_) {
      // update the trace formula with the consecutive axioms
      // needed for it to be unsat with all the nonconsecutive axioms
      // it will be updated again later anyway
      for (auto ax : consecutive_axioms) {
        size_t max_k = abs_ts_.only_curr(ax) ? reached_k_ + 1
                                                : reached_k_;
        for (size_t k = 0; k <= max_k; ++k) {
          abs_bmc_formula = super::solver_->make_term(
              And, abs_bmc_formula, abs_unroller_.at_time(ax, k));
        }
      }

      nonconsecutive_axioms =
          reduce_nonconsecutive_axioms(abs_bmc_formula, nonconsecutive_axioms);
    }

    // First collect all the indices used in nonconsecutive axioms
    // these will need to be prophecized
    UnorderedTermSet instantiations;
    for (AxiomInstantiation ax_inst : nonconsecutive_axioms) {
      assert(ax_inst.instantiations.size()
             == 1);  // expecting only one index instantiated
      for (auto inst : ax_inst.instantiations) {
        instantiations.insert(inst);
      }
    }

    // vector of pairs
    // first: prophecy variable
    // second: target (a history variable for non-zero delay)
    vector<pair<Term, Term>> proph_vars;
    for (auto timed_idx : instantiations) {
      // number of steps before the property violation
      size_t delay =
          reached_k_ + 1 - abs_unroller_.get_curr_time(timed_idx);
      // Prophecy Modifier will add prophecy and history variables
      // automatically here but it does NOT update the property
      Term idx = abs_unroller_.untime(timed_idx);
      // can't target a non-current state variable
      // because the target will appear in the updated property
      assert(delay > 0 || abs_ts_.only_curr(idx));
      proph_vars.push_back(pm_.get_proph(idx, delay));
    }

    assert(instantiations.size() == proph_vars.size());
    logger.log(1, "Added {} prophecy variables", proph_vars.size());

    // now update bad_ and add the prophecy variables to the index set
    // the property would be updated as
    // (proph1=target1 /\ ... /\ prophn=targetn) -> prop
    // but because we're working with bad, this is equivalent to
    // (proph1=target1 /\ ... /\ prophn=targetn) /\ bad
    // where bad = !prop
    for (auto p : proph_vars) {
      Term proph_var = p.first;
      Term target = p.second;
      aae_.add_index(proph_var);
      super::bad_ = super::solver_->make_term(
          And,
          super::solver_->make_term(Equal, proph_var, target),
          super::bad_);
    }

    // need to update the bmc formula with the transformations
    // to abs_ts_
    abs_bmc_formula = get_bmc_formula(reached_k_ + 1);

    // search for axioms again but don't include nonconsecutive ones
    bool ok =
        aae_.enumerate_axioms(abs_bmc_formula, reached_k_ + 1, false);
    // should be guaranteed to rule out counterexamples at this bound
    assert(ok);
    consecutive_axioms = aae_.get_consecutive_axioms();
    assert(!aae_.get_nonconsecutive_axioms().size());
  }

  if (super::options_.cegp_axiom_red_ && consecutive_axioms.size()) {
    reduce_consecutive_axioms(abs_bmc_formula, consecutive_axioms);
  }

  if (consecutive_axioms.size() > 0) {
    refine_ts(consecutive_axioms);
    num_added_axioms_ += consecutive_axioms.size();
    logger.log(1, "CEGP: refine added {} axiom(s)", num_added_axioms_);
  }

  // able to successfully refine
  return true;
}

// helpers
template <class Prover_T>
Term CegProphecyArrays<Prover_T>::get_bmc_formula(size_t b)
{
  Term abs_bmc_formula = abs_unroller_.at_time(abs_ts_.init(), 0);
  for (int k = 0; k < b; ++k) {
    abs_bmc_formula = super::solver_->make_term(
        And, abs_bmc_formula, abs_unroller_.at_time(abs_ts_.trans(), k));
  }

  return super::solver_->make_term(
      And, abs_bmc_formula, abs_unroller_.at_time(super::bad_, b));
}

template <class Prover_T>
void CegProphecyArrays<Prover_T>::reduce_consecutive_axioms(
    const Term & abs_bmc_formula, UnorderedTermSet & consec_ax)
{
  UnorderedTermSet reduced_ax;
  super::solver_->push();

  super::solver_->assert_formula(abs_bmc_formula);

  TermVec assumps;
  UnorderedTermMap label2ax;
  Term unrolled_ax;
  Term lbl;
  for (auto ax : consec_ax) {
    unrolled_ax = super::solver_->make_term(true);
    size_t max_k =
        abs_ts_.only_curr(ax) ? reached_k_ + 1 : reached_k_;
    for (size_t k = 0; k <= max_k; ++k) {
      unrolled_ax = super::solver_->make_term(
          And, unrolled_ax, abs_unroller_.at_time(ax, k));
    }

    lbl = label(unrolled_ax);
    label2ax[lbl] = ax;
    super::solver_->assert_formula(
        super::solver_->make_term(Implies, lbl, unrolled_ax));
    assumps.push_back(lbl);
  }

  Result res = super::solver_->check_sat_assuming(assumps);
  assert(res.is_unsat());

  UnorderedTermSet core;
  super::solver_->get_unsat_core(core);

  for (auto l : assumps) {
    if (core.find(l) == core.end()) {
      // if not in core, then remove from axioms
      size_t num_erased = consec_ax.erase(label2ax.at(l));
      assert(num_erased);  // expecting axiom to be in set
    }
  }

  super::solver_->pop();
}

template <class Prover_T>
AxiomVec CegProphecyArrays<Prover_T>::reduce_nonconsecutive_axioms(
    const Term & abs_bmc_formula, const AxiomVec & nonconsec_ax)
{
  // map from delay to the target (over ts vars) and a vector of axioms using
  // that target using to sort: rely on sortedness of map to put in ascending
  // order of delay
  map<int, unordered_map<Term, AxiomVec>> map_nonconsec_ax;
  Term unrolled_idx;
  Term idx;
  for (auto ax_inst : nonconsec_ax) {
    // expecting only a single index to be instantiated
    assert(ax_inst.instantiations.size() == 1);
    unrolled_idx = *(ax_inst.instantiations.begin());
    size_t delay =
        reached_k_ + 1 - abs_unroller_.get_curr_time(unrolled_idx);
    idx = abs_unroller_.untime(unrolled_idx);
    map_nonconsec_ax[delay][idx].push_back(ax_inst);
  }

  vector<AxiomVec> sorted_nonconsec_ax;
  for (auto elem1 : map_nonconsec_ax) {
    for (auto elem2 : elem1.second) {
      sorted_nonconsec_ax.push_back(elem2.second);
    }
  }

  // now try to drop the targets with larger delay one-by-one
  AxiomVec red_nonconsec_ax;
  // a vector of bools indicating if each AxiomVec should be included
  // starts with including all
  vector<bool> include_axioms(sorted_nonconsec_ax.size(), true);
  super::solver_->push();

  Result res;
  for (int i = sorted_nonconsec_ax.size() - 1; i >= 0; --i) {
    super::solver_->push();
    // assert all the included axioms except the i-th
    for (size_t j = 0; j < include_axioms.size(); j++) {
      if (j == i || !include_axioms[j]) {
        continue;
      }
      for (auto ax_inst : sorted_nonconsec_ax[j]) {
        super::solver_->assert_formula(ax_inst.ax);
      }
    }

    res = super::solver_->check_sat();
    if (res.is_unsat()) {
      // the i-th vector of axioms is unneeded
      include_axioms[i] = false;
    }
    super::solver_->pop();
  }

  super::solver_->pop();

  assert(include_axioms.size() == sorted_nonconsec_ax.size());
  for (size_t i = 0; i < include_axioms.size(); ++i) {
    if (include_axioms[i]) {
      red_nonconsec_ax.insert(red_nonconsec_ax.end(),
                              sorted_nonconsec_ax[i].begin(),
                              sorted_nonconsec_ax[i].end());
    }
  }

  return red_nonconsec_ax;
}

template <class Prover_T>
Term CegProphecyArrays<Prover_T>::label(const Term & t)
{
  auto it = labels_.find(t);
  if (it != labels_.end()) {
    return labels_.at(t);
  }

  unsigned i = 0;
  Term l;
  while (true) {
    try {
      l = super::solver_->make_symbol(
          "assump_" + std::to_string(t->hash()) + "_" + std::to_string(i),
          super::solver_->make_sort(BOOL));
      break;
    }
    catch (IncorrectUsageException & e) {
      ++i;
    }
    catch (SmtException & e) {
      throw e;
    }
  }

  labels_[t] = l;
  return l;
}

template <class Prover_T>
void CegProphecyArrays<Prover_T>::refine_ts(const UnorderedTermSet & consecutive_axioms)
{
  // add consecutive axioms to the system
  // TODO: make sure we're adding current / next correctly
  RelationalTransitionSystem & rts =
    static_cast<RelationalTransitionSystem &>(abs_ts_);
  for (const auto & ax : consecutive_axioms) {
    if (reached_k_ == -1) {
      // if only checking initial state
      // need to add to init
      rts.constrain_init(ax);
    }

    rts.constrain_trans(ax);
    if (rts.only_curr(ax)) {
      // add the next state version if it's an invariant over current state vars
      rts.constrain_trans(abs_ts_.next(ax));
    }
  }

  refine_subprover_ts(consecutive_axioms);
}

template <class Prover_T>
void CegProphecyArrays<Prover_T>::refine_subprover_ts(const UnorderedTermSet & consecutive_axioms)
{
  // No-Op
}

template <>
void CegProphecyArrays<IC3IA>::refine_subprover_ts(const UnorderedTermSet & consecutive_axioms)
{
  const RelationalTransitionSystem & rts =
    static_cast<const RelationalTransitionSystem &>(abs_ts_);
  RelationalTransitionSystem & sub_rts =
    static_cast<RelationalTransitionSystem &>(ts_);

  // add predicates from init and bad
  UnorderedTermSet preds;
  get_predicates(solver_, abs_ts_.init(), preds, false, false, true);
  get_predicates(solver_, bad_, preds, false, false, true);
  // instead of add previously found predicates, we add all the predicates in frame 1
  // preds.insert(predset_.begin(), predset_.end());
  get_predicates(solver_, get_frame_term(1), preds, false, false, true);
  predset_.clear();
  predlbls_.clear();

  // don't add boolean symbols that are never used in the system
  // this is an optimization and a fix for some options
  // if using mathsat with bool_model_generation
  // it will fail to get the value of symbols that don't
  // appear in the query
  // thus we don't include those symbols in our cubes
  UnorderedTermSet used_symbols;
  get_free_symbolic_consts(ts_.init(), used_symbols);
  get_free_symbolic_consts(ts_.trans(), used_symbols);
  get_free_symbolic_consts(bad_, used_symbols);

  // reset init and trans -- done with calling ia_.do_abstraction
  // then add all boolean constants as (precise) predicates
  for (const auto & p : ia_.do_abstraction()) {
    assert(p->is_symbolic_const());
    if (used_symbols.find(p) != used_symbols.end()) {
      preds.insert(p);
    }
  }

  // reset the solver
  reset_solver();

  // add predicates
  for (const auto &p : preds) {
    add_predicate(p);
  }
}

// ceg-prophecy is incremental for ic3ia
template class CegProphecyArrays<IC3IA>;

// the below engines work when ceg-prophecy is non-incremental.
template class CegProphecyArrays<Bmc>;
template class CegProphecyArrays<BmcSimplePath>;
template class CegProphecyArrays<KInduction>;
template class CegProphecyArrays<InterpolantMC>;
template class CegProphecyArrays<ModelBasedIC3>;

#ifdef WITH_MSAT_IC3IA
template class CegProphecyArrays<MsatIC3IA>;
#endif
}  // namespace pono
