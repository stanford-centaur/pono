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

#include <cassert>

#include "modifiers/static_coi.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/term_walkers.h"

using namespace smt;

namespace pono {

OpInpAbstractor::OpInpAbstractor(
    const TransitionSystem & conc_ts,
    TransitionSystem & abs_ts,
    const OpSet & op_to_abstract,
    const smt::Term & prop,  // it is okay to use bad
    int verbosity)
    : OpAbstractor(conc_ts, abs_ts)
{
  abstract_ts(conc_ts, abs_ts, op_to_abstract, prop, verbosity);
  unroller_ = std::make_unique<Unroller>(abs_ts);
}  // OpAbstractor

void OpInpAbstractor::abstract_ts(
    const TransitionSystem & in_ts,
    TransitionSystem & out_ts,
    const OpSet & op_to_abstract,
    const smt::Term & prop,  // it is okay to use bad
    int verbosity)
{
  // copy if not done outside
  if (out_ts.inputvars().empty() && out_ts.statevars().empty()) out_ts = in_ts;

  TermOpCollector op_collector(out_ts.solver());

  UnorderedTermSet term_op_out;
  for (const auto & s_update : out_ts.state_updates()) {
    op_collector.find_matching_terms(
        s_update.second, op_to_abstract, term_op_out);
  }  // walk state update functions

  unsigned dummy_input_cnt = 0;
  UnorderedTermMap replacement;
  for (const auto & t : term_op_out) {
    // add to <op, result, arg>
    op_abstracted.push_back(OpAbstract());
    op_abstracted.back().op = t->get_op();
    op_abstracted.back().args = TermVec(t->begin(), t->end());
    op_abstracted.back().original = t;

    bool succ = false;
    do {
      try {
        op_abstracted.back().result = out_ts.make_inputvar(
            "_dummy_input_cnt_" + std::to_string(dummy_input_cnt++),
            t->get_sort());
        succ = true;
        dummy_inputs_.insert(op_abstracted.back().result);
      }
      catch (PonoException & e) {
      }
    } while (!succ);

    replacement.emplace(t, op_abstracted.back().result);
  }  // for each term

  out_ts.replace_terms(replacement);
  {
    StaticConeOfInfluence coi(out_ts, { prop }, verbosity);
  }
}

bool check_possible(const Term & res,
                    const TermVec & arg,
                    Op op,
                    const SmtSolver & orig_solver)
{
  SmtSolver s = create_solver(BTOR, false, false, false);
  TermTranslator tr(s);
  TermVec tr_arg;

  assert(res->is_value());
  for (const auto & val : arg) {
    assert(val->is_value());
    tr_arg.push_back(tr.transfer_term(val));
  }
  auto rhs = s->make_term(op, tr_arg);
  s->assert_formula(s->make_term(Equal, rhs, tr.transfer_term(res)));
  auto eq_sat = s->check_sat();
  return eq_sat.is_sat();
}

bool OpInpAbstractor::refine_with_constraints(
    // TODO: how to put the trace in here?
    const smt::TermVec & cexs,
    const Term & bad,
    smt::TermVec & out)
{
  // short-cut
  if (op_abstracted.empty()) return false;

  assert(cexs.size() > 1);

  // unroll from 0 --> ...
  size_t trace_pos = 0;

  const auto & solver_ = abs_ts_.solver();
  solver_->push();
  bool refined = false;

  while (trace_pos != cexs.size()) {
    refined = false;
    const auto & curr_term = cexs.at(trace_pos);
    logger.log(3, "Refine cex t{} : {}", trace_pos, curr_term->to_string());

    if (trace_pos == 0) {
      solver_->assert_formula(unroller_->at_time(abs_ts_.init(), 0));
      solver_->assert_formula(unroller_->at_time(curr_term, 0));
    } else {
      // trans will include constraints
      solver_->assert_formula(
          unroller_->at_time(abs_ts_.trans(), trace_pos - 1));
      solver_->assert_formula(unroller_->at_time(curr_term, trace_pos));
      auto res = solver_->check_sat();
      assert(res.is_sat());

      // find term (result) in the model
      UnorderedTermSet vars;
      const auto & prev_term = cexs.at(trace_pos - 1);
      get_free_symbolic_consts(prev_term, vars);
      for (auto & abs_term : op_abstracted) {
        if (vars.find(abs_term.result) != vars.end()) {
          // now eval the result (trace_pos-1) and the args (trace_pos-1)
          // and see if violated
          auto res_val = solver_->get_value(
              unroller_->at_time(abs_term.result, trace_pos - 1));
          TermVec arg_val;
          for (const auto & arg : abs_term.args) {
            arg_val.push_back(
                solver_->get_value(unroller_->at_time(arg, trace_pos - 1)));
          }
          if (!check_possible(res_val, arg_val, abs_term.op, solver_)) {
            // create an assumption
            Term expected_res_val = solver_->make_term(abs_term.op, arg_val);
            logger.log(3,
                       "Refine {} ({} {} {}) = {} =/= {}",
                       abs_term.op.to_string(),
                       arg_val.size() > 0 ? arg_val.at(0)->to_string() : "",
                       arg_val.size() > 1 ? arg_val.at(1)->to_string() : "",
                       arg_val.size() > 2 ? "..." : "",
                       expected_res_val->to_string(),
                       res_val->to_string());
            auto refine_res =
                solver_->make_term(Equal, abs_term.result, expected_res_val);

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

            abs_term.refine_count++;
            refined = true;
            // you may want to break here
          } else {
            logger.log(3,
                       "NO Refine {} ({} {} {}) = {}",
                       abs_term.op.to_string(),
                       arg_val.size() > 0 ? arg_val.at(0)->to_string() : "",
                       arg_val.size() > 1 ? arg_val.at(1)->to_string() : "",
                       arg_val.size() > 2 ? "..." : "",
                       res_val->to_string());
          }
          // check if this is possible
        }  // end if find this abs_term exists
      }
    }  // else (trace_pos > 0)

    if (refined) break;

    ++trace_pos;
  }  // end of while

  assert(trace_pos > 0);

  solver_->pop();
  return refined;
}  // end of OpAbstractor::refine_with_constraints

// ----------------------------------------------------------------

OpUfAbstractor::OpUfAbstractor(const TransitionSystem & conc_ts,
                               TransitionSystem & abs_ts,
                               const OpSet & op_to_abstract)
    : OpAbstractor(conc_ts, abs_ts), uf_extractor_(abs_ts.solver())
{
  abstract_ts(conc_ts, abs_ts, op_to_abstract);
  unroller_ = std::make_unique<Unroller>(abs_ts);
}  // OpAbstractor

void OpUfAbstractor::abstract_ts(const TransitionSystem & in_ts,
                                 TransitionSystem & out_ts,
                                 const OpSet & op_to_abstract)
{
  // copy if not done outside
  if (out_ts.inputvars().empty() && out_ts.statevars().empty()) out_ts = in_ts;
  const auto & solver_ = out_ts.solver();

  TermOpCollector op_collector(out_ts.solver());

  UnorderedTermSet term_op_out;
  for (const auto & s_update : out_ts.state_updates()) {
    op_collector.find_matching_terms(
        s_update.second, op_to_abstract, term_op_out);
  }  // walk state update functions

  unsigned dummy_uf_cnt = 0;
  UnorderedTermMap replacement;
  for (const auto & t : term_op_out) {
    // compute a signature of function
    op_abstracted.push_back(OpUfAbstract());
    auto & abitem = op_abstracted.back();
    abitem.op = t->get_op();
    abitem.args = TermVec(t->begin(), t->end());
    abitem.original = t;

    std::string func_signature = t->get_op().to_string();
    for (const auto & a : abitem.args) {
      func_signature += "_" + a->get_sort()->to_string();
    }

    auto pos = uf_set_.find(func_signature);
    if (pos == uf_set_.end()) {
      SortVec func_sort;
      for (const auto & a : abitem.args) func_sort.push_back(a->get_sort());
      func_sort.push_back(t->get_sort());  // the last one is the result sort

      // now make function
      Sort funsort = solver_->make_sort(FUNCTION, func_sort);

      bool new_uf_succ = false;
      Term f;
      while (!new_uf_succ) {
        try {
          std::string new_name = "__uf" + std::to_string(dummy_uf_cnt++);
          f = solver_->make_symbol(new_name, funsort);
          new_uf_succ = true;
        }
        catch (...) {
        }
      }
      assert(f != nullptr);

      pos = uf_set_.emplace(func_signature, f).first;
    }  // end if not found
    abitem.uf = pos->second;

    TermVec arg;  // f arg1 arg2 ...
    arg.push_back(pos->second);
    arg.insert(arg.end(), abitem.args.begin(), abitem.args.end());

    abitem.result = out_ts.make_term(Apply, arg);
    replacement.emplace(t, abitem.result);
  }  // end for each found op term

  out_ts.replace_terms(replacement);
}  // end of OpUfAbstractor::abstract_ts

bool OpUfAbstractor::refine_with_constraints(
    // TODO: how to put the trace in here?
    const smt::TermVec & cexs,
    const Term & bad,
    smt::TermVec & out)
{
  // short-cut
  if (op_abstracted.empty()) return false;

  assert(cexs.size() > 1);

  // unroll from 0 --> ...
  size_t trace_pos = 0;

  const auto & solver_ = abs_ts_.solver();
  solver_->push();
  bool refined = false;

  while (trace_pos != cexs.size()) {
    refined = false;
    const auto & curr_term = cexs.at(trace_pos);
    logger.log(3, "Refine cex t{} : {}", trace_pos, curr_term->to_string());

    if (trace_pos == 0) {
      solver_->assert_formula(unroller_->at_time(abs_ts_.init(), 0));
      solver_->assert_formula(unroller_->at_time(curr_term, 0));
    } else {
      // trans will include constraints
      solver_->assert_formula(
          unroller_->at_time(abs_ts_.trans(), trace_pos - 1));
      solver_->assert_formula(unroller_->at_time(curr_term, trace_pos));
      auto res = solver_->check_sat();
      assert(res.is_sat());

      // find term (result) in the model
      UnorderedTermSet ufs;
      const auto & prev_term = cexs.at(trace_pos - 1);
      // get_free_symbols( prev_ptr->target.term, vars);
      // here we need ptr->target.term 's prev expression
      // extract all f,  find the existence of some f?
      uf_extractor_.extract_uf_in_term(
          solver_->substitute(prev_term, abs_ts_.state_updates()), ufs);

      for (auto & abs_term : op_abstracted) {
        if (ufs.find(abs_term.uf) != ufs.end()) {
          // now eval the result (trace_pos-1) and the args (trace_pos-1)
          // and see if violated
          auto res_val = solver_->get_value(
              unroller_->at_time(abs_term.result, trace_pos - 1));
          TermVec arg_val;
          for (const auto & arg : abs_term.args) {
            arg_val.push_back(
                solver_->get_value(unroller_->at_time(arg, trace_pos - 1)));
          }
          if (!check_possible(res_val, arg_val, abs_term.op, solver_)) {
            // create an assumption
            Term expected_res_val = solver_->make_term(abs_term.op, arg_val);
            TermVec apply_farg;
            apply_farg.push_back(abs_term.uf);
            apply_farg.insert(apply_farg.end(), arg_val.begin(), arg_val.end());

            Term uf_expr = solver_->make_term(Apply, apply_farg);
            logger.log(3,
                       "Refine {} ({} {} {}) = {} =/= {}",
                       abs_term.op.to_string(),
                       arg_val.size() > 0 ? arg_val.at(0)->to_string() : "",
                       arg_val.size() > 1 ? arg_val.at(1)->to_string() : "",
                       arg_val.size() > 2 ? "..." : "",
                       expected_res_val->to_string(),
                       res_val->to_string());
            auto refine_res =
                solver_->make_term(Equal, uf_expr, expected_res_val);
            out.push_back(refine_res);

            abs_term.refine_count++;
            refined = true;
            break;
            // you may want to break here
          } else {
            logger.log(3,
                       "NO Refine {} ({} {} {}) = {}",
                       abs_term.op.to_string(),
                       arg_val.size() > 0 ? arg_val.at(0)->to_string() : "",
                       arg_val.size() > 1 ? arg_val.at(1)->to_string() : "",
                       arg_val.size() > 2 ? "..." : "",
                       res_val->to_string());
          }
          // check if this is possible
        }  // end if find this abs_term exists
      }
    }  // else (trace_pos > 0)

    if (refined) break;

    ++trace_pos;
  }  // end of while

  assert(trace_pos > 0);

  solver_->pop();
  return refined;
}  // end of OpUfAbstractor::refine_with_constraints

// ------------------------------------------

void UfExtractor::extract_uf_in_term(const smt::Term & t,
                                     smt::UnorderedTermSet & uf_out)
{
  out_ = &uf_out;
  Term tmp(t);
  visit(tmp);
  out_ = NULL;
}

WalkerStepResult UfExtractor::visit_term(smt::Term & term)
{
  assert(out_);

  if (preorder_) {
    Op op = term->get_op();
    if (op.prim_op == Apply) {
      assert(term->begin() != term->end());
      out_->insert(*(term->begin()));
    }
  }

  return Walker_Continue;
}

}  // namespace pono
