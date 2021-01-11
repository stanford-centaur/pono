/*********************                                                  */
/*! \file op_abstractor.cpp
** \verbatim
** Top contributors (to current version):
**   Hongce Zhang
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Abstract complex operators with uninterpreted functions
**
**
**
**/

#include "modifiers/op_abstractor.h"
#include "modifiers/static_coi.h"
#include "utils/term_walkers.h"
#include "core/unroller.h"
#include "smt/available_solvers.h"

#include <cassert>

using namespace smt;

namespace pono {

OpAbstractor::OpAbstractor(
    const TransitionSystem & conc_ts,
    TransitionSystem & abs_ts,
    const OpSet & op_to_abstract,
    const smt::Term & prop, // it is okay to use bad
    int verbosity
    ) : Abstractor(conc_ts, abs_ts)
{ // copy if not done outside
  if (abs_ts.inputvars().empty() && abs_ts.statevars().empty())
    abs_ts = conc_ts;

  TermOpCollector op_collector(abs_ts.solver());

  UnorderedTermSet term_op_out;
  for (const auto & s_update : abs_ts.state_updates()) {
    op_collector.find_matching_terms(s_update.second, op_to_abstract, term_op_out);
  } // walk state update functions

  unsigned dummy_input_cnt = 0;
  UnorderedTermMap replacement;
  for (const auto & t : term_op_out) {
    // add to <op, result, arg>
    op_abstracted.push_back(OpAbstract());
    op_abstracted.back().op = t->get_op();
    op_abstracted.back().args = TermVec(t->begin(), t->end());
    op_abstracted.back().original = t;

    bool succ = false;
    do{
      try {
        op_abstracted.back().result =
          abs_ts.make_inputvar(
            "_dummy_input_cnt_" + std::to_string(dummy_input_cnt++),
            t->get_sort());
        succ = true;
      } catch (PonoException & e) {
      }
    } while(!succ);

    replacement.emplace(t, op_abstracted.back().result);
  } // for each term

  abs_ts.replace_terms(replacement);
  {
    StaticConeOfInfluence coi(abs_ts, { prop }, verbosity);
  } 
} // OpAbstractor

bool check_possible(const Term & res, const TermVec & arg, Op op, Term & rhs_val,
    const SmtSolver & orig_solver) {
  SmtSolver s = create_solver(BTOR, false, false ,false);
  TermTranslator tr(s);
  TermTranslator tr_back(orig_solver);
  TermVec tr_arg;

  assert(res->is_value());
  for (const auto & val : arg) {
    assert(val->is_value());
    tr_arg.push_back(tr.transfer_term(val));
  }
  auto rhs = s->make_term(op, tr_arg);
  s->assert_formula(
    s->make_term(
      Equal,
      rhs, 
      tr.transfer_term(res))
  );
  auto eq_sat = s->check_sat();
  if (eq_sat.is_unsat()) {
    rhs_val = tr_back.transfer_term(s->get_value(rhs));
  }
  return eq_sat.is_sat();
}

bool OpAbstractor::refine_with_constraints(
    // TODO: how to put the trace in here?
    const ProofGoal * goal_at_init,
    smt::TermVec & out) {
  
  assert(goal_at_init);
  assert(goal_at_init->idx == 0);

  Unroller unroller(abs_ts_);
  // unroll from 0 --> ...
  const ProofGoal * ptr = goal_at_init;
  const ProofGoal * prev_ptr = NULL;
  unsigned time = 0;

  const auto & solver_ = abs_ts_.solver();
  solver_->push();
  while(ptr) {
    bool refined = false;
    if (time == 0) {
      solver_->assert_formula( unroller.at_time(abs_ts_.init(), 0) );
      solver_->assert_formula( unroller.at_time(ptr->target.term, 0));
    } else {
      // trans will include constraints
      solver_->assert_formula( unroller.at_time(abs_ts_.trans(), time-1) );
      solver_->assert_formula( unroller.at_time(ptr->target.term, time));
      auto res = solver_->check_sat();
      assert (res.is_sat());

      // find term (result) in the model
      UnorderedTermSet vars;
      assert(prev_ptr);
      get_free_symbols( prev_ptr->target.term, vars);
      for (auto & abs_term : op_abstracted) {
        if (vars.find(abs_term.result) != vars.end()) {
          // now eval the result (time-1) and the args (time-1)
          // and see if violated
          auto res_val = solver_->get_value(
            unroller.at_time(abs_term.result, time-1));
          TermVec arg_val;
          for (const auto & arg : abs_term.args) {
            arg_val.push_back(
              solver_->get_value(
                unroller.at_time(arg, time-1)));
          }
          Term res_val;
          if (!check_possible(res_val, arg_val, abs_term.op, res_val, solver_)) {
            // create an assumption
            auto refine_res = solver_->make_term(Equal, abs_term.result , res_val);

            unsigned arg_idx = 0;
            Term func_arg;
            for (const auto & arg : abs_term.args) {
              auto eq = solver_->make_term(Equal, arg, arg_val.at(arg_idx++));
              if (func_arg == nullptr)
                func_arg = eq;
              else
                func_arg = solver_->make_term(And, func_arg, eq);
            }
            auto cnstr = solver_->make_term(Implies, func_arg, refine_res);
            out.push_back(cnstr);
            
            abs_term.refine_count ++;
            refined = true;
            // you may want to break here
          }
          // check if this is possible
        } // end if find this abs_term exists
      }
    } // else (time > 0)

    if (refined)
      break;

    ++ time;
    prev_ptr = ptr;
    ptr = ptr->next;
  } // end of while

  assert(time > 0);

  solver_->pop();
} // end of OpAbstractor::refine_with_constraints


} // namespace pono
