/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann
 ** This file is part of the pono project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief
 **
 **
 **/

#include "btor2_encoder.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>

#include "smt-switch/utils.h"
#include "utils/logger.h"

using namespace smt;
using namespace std;

namespace pono {

// Maps for use in conversion
const unordered_map<Btor2Tag, int> basemap({ { BTOR2_TAG_const, 2 },
                                             { BTOR2_TAG_constd, 10 },
                                             { BTOR2_TAG_consth, 16 } });

const unordered_map<Btor2Tag, smt::PrimOp> bvopmap({
    { BTOR2_TAG_add, BVAdd },
    { BTOR2_TAG_and, BVAnd },
    // { BTOR2_TAG_bad, },
    { BTOR2_TAG_concat, Concat },
    //{ BTOR2_TAG_const, },
    //{ BTOR2_TAG_constraint, },
    //{ BTOR2_TAG_constd, },
    //{ BTOR2_TAG_consth, },
    //{ BTOR2_TAG_dec, },
    // { BTOR2_TAG_eq, BVComp }, // handled this specially, because could also
    // have array arguments
    //{ BTOR2_TAG_fair, },
    { BTOR2_TAG_iff, BVComp },
    // { BTOR2_TAG_implies, Implies }, // handle specially (needs to work with
    // boolector AND other solvers), gets complicated with bools and BV(1) are
    // aliased
    //{ BTOR2_TAG_inc, },
    //{ BTOR2_TAG_init, },
    //{ BTOR2_TAG_input, },
    //{ BTOR2_TAG_ite, Ite },
    //{ BTOR2_TAG_justice, },
    { BTOR2_TAG_mul, BVMul },
    { BTOR2_TAG_nand, BVNand },
    // handled this specially, because could also have array arguments
    // { BTOR2_TAG_neq, Distinct },
    { BTOR2_TAG_neg, BVNeg },
    //{ BTOR2_TAG_next, },
    { BTOR2_TAG_nor, BVNor },
    { BTOR2_TAG_not, BVNot },
    //{ BTOR2_TAG_one, },
    //{ BTOR2_TAG_ones, },
    { BTOR2_TAG_or, BVOr },
    //{ BTOR2_TAG_output, },
    // { BTOR2_TAG_read, Select }, // handle specially -- make sure it's casted
    // to bv
    //{ BTOR2_TAG_redand, },
    //{ BTOR2_TAG_redor, },
    //{ BTOR2_TAG_redxor, },
    // { BTOR2_TAG_rol, },
    // { BTOR2_TAG_ror, },
    //{ BTOR2_TAG_saddo, },
    { BTOR2_TAG_sdiv, BVSdiv },
    //{ BTOR2_TAG_sdivo, },
    //{ BTOR2_TAG_sext, },
    { BTOR2_TAG_sgt, BVSgt },
    { BTOR2_TAG_sgte, BVSge },
    //{ BTOR2_TAG_slice, },
    { BTOR2_TAG_sll, BVShl },
    { BTOR2_TAG_slt, BVSlt },
    { BTOR2_TAG_slte, BVSle },
    //{ BTOR2_TAG_sort, },
    { BTOR2_TAG_smod, BVSmod },
    //{ BTOR2_TAG_smulo, },
    { BTOR2_TAG_sra, BVAshr },
    { BTOR2_TAG_srem, BVSrem },
    { BTOR2_TAG_srl, BVLshr },
    //{ BTOR2_TAG_ssubo, },
    //{ BTOR2_TAG_state, },
    { BTOR2_TAG_sub, BVSub },
    //{ BTOR2_TAG_uaddo, },
    { BTOR2_TAG_udiv, BVUdiv },
    //{ BTOR2_TAG_uext, },
    { BTOR2_TAG_ugt, BVUgt },
    { BTOR2_TAG_ugte, BVUge },
    { BTOR2_TAG_ult, BVUlt },
    { BTOR2_TAG_ulte, BVUle },
    //{ BTOR2_TAG_umulo, },
    { BTOR2_TAG_urem, BVUrem },
    //{ BTOR2_TAG_usubo, },
    //{ BTOR2_TAG_write, Store }, // handle specially -- make sure it's casted
    // to bv
    { BTOR2_TAG_xnor, BVXnor },
    { BTOR2_TAG_xor, BVXor },
    //{ BTOR2_TAG_zero, }
});

const unordered_map<Btor2Tag, smt::PrimOp> boolopmap({
    { BTOR2_TAG_and, And },
    { BTOR2_TAG_or, Or },
    { BTOR2_TAG_xor, Xor },
    { BTOR2_TAG_not, Not },
    // { BTOR2_TAG_implies, Implies },
    { BTOR2_TAG_iff, Equal },
    // handling specially -- could have array arguments
    // { BTOR2_TAG_eq, Equal },
    //{ BTOR2_TAG_neq, Distinct }
});

Term BTOR2Encoder::bool_to_bv(const Term & t) const
{
  if (t->get_sort()->get_sort_kind() == BOOL) {
    Sort bv1sort = solver_->make_sort(BV, 1);
    return solver_->make_term(
        Ite, t, solver_->make_term(1, bv1sort), solver_->make_term(0, bv1sort));
  } else {
    return t;
  }
}

Term BTOR2Encoder::bv_to_bool(const Term & t) const
{
  Sort sort = t->get_sort();
  if (sort->get_sort_kind() == BV) {
    if (sort->get_width() != 1) {
      throw PonoException("Can't convert non-width 1 bitvector to bool.");
    }
    return solver_->make_term(
        Equal, t, solver_->make_term(1, solver_->make_sort(BV, 1)));
  } else {
    return t;
  }
}

TermVec BTOR2Encoder::lazy_convert(const TermVec & tvec) const
{
  TermVec res;
  res.reserve(tvec.size());

  unsigned int num_bools = 0;
  bool wide_bvs;
  Sort sort;
  UnorderedSortSet sortset;
  for (auto t : tvec) {
    res.push_back(t);

    sort = t->get_sort();
    sortset.insert(sort);

    if (sort->get_sort_kind() == BOOL) {
      num_bools++;
    } else if (!(sort->get_sort_kind() == BV && sort->get_width() == 1)) {
      wide_bvs = true;
    }
  }

  if (sortset.size() > 1) {
    if (num_bools > tvec.size() / 2 && !wide_bvs) {
      for (size_t i = 0; i < res.size(); i++) {
        res[i] = bv_to_bool(res[i]);
      }
    } else {
      for (size_t i = 0; i < res.size(); i++) {
        res[i] = bool_to_bv(res[i]);
      }
    }
  }

  return res;
}

// to handle the case where yosys generate sth. like this
// state (with no name)
// output (with name)
// for Verilog :  output reg xxx;

// this function go over the file and record this case
// when we encounter it again we shall replace the state's name
void BTOR2Encoder::preprocess(const std::string & filename)
{
  FILE * input_file = fopen(filename.c_str(), "r");

  if (!input_file) {
    throw PonoException("Could not open " + filename);
  }

  Btor2Parser * reader = btor2parser_new();

  if (!btor2parser_read_lines(reader, input_file)) {
    throw PonoException(std::string(btor2parser_error(reader)));
  }

  Btor2LineIterator bt2_it = btor2parser_iter_init(reader);

  std::unordered_set<uint64_t> unamed_state_ids;
  while (Btor2Line * bt2_line = btor2parser_iter_next(&bt2_it)) {
    if (bt2_line->tag == BTOR2_TAG_state) {
      if (!bt2_line->symbol) {  // if we see state has no name, record it
        unamed_state_ids.insert(bt2_line->id);
      }
    } else if (bt2_line->tag == BTOR2_TAG_output && bt2_line->symbol) {
      // if we see an output with name
      if (bt2_line->nargs == 0)
        throw PonoException("Missing term for id "
                            + std::to_string(bt2_line->id));
      // we'd like to know if it refers to a state without name
      auto pos =
          unamed_state_ids.find(bt2_line->args[0]);  // so, *pos is its btor id
      if (pos != unamed_state_ids.end()) {
        // in such case, we record that we can name it with this new name
        auto state_name_pos = state_renaming_table_.find(*pos);
        if (state_name_pos == state_renaming_table_.end()) {
          state_renaming_table_.insert(std::make_pair(*pos, bt2_line->symbol));
        }  // otherwise we already have a name for it, then just ignore this one
      }
    }  // end of if input
  }  // end of while

  fclose(input_file);
  btor2parser_delete(reader);
}  // end of preprocess

void BTOR2Encoder::parse(const std::string filename)
{
  FILE * input_file = fopen(filename.c_str(), "r");

  if (!input_file) {
    throw PonoException("Could not open " + filename);
  }

  Btor2Parser * reader = btor2parser_new();

  if (!btor2parser_read_lines(reader, input_file)) {
    throw PonoException(std::string(btor2parser_error(reader)));
  }

  // local variables
  uint64_t num_states = 0;
  std::unordered_map<int64_t, uint64_t> id2statenum;
  smt::Sort linesort;
  smt::TermVec termargs;
  std::string orig_symbol;
  std::string new_symbol;

  Btor2LineIterator bt2_it = btor2parser_iter_init(reader);
  while (Btor2Line * bt2_line = btor2parser_iter_next(&bt2_it)) {
    /******************************** Identify sort
     * ********************************/
    if (bt2_line->tag != BTOR2_TAG_sort && bt2_line->sort.id) {
      linesort = sorts_.at(bt2_line->sort.id);
    }

    /******************************** Gather term arguments
     * ********************************/
    termargs.clear();
    termargs.reserve(bt2_line->nargs);
    for (size_t i_arg = 0; i_arg < bt2_line->nargs; i_arg++) {
      bool negated = false;
      int64_t arg_bt2_id = bt2_line->args[i_arg];
      if (arg_bt2_id < 0) {
        negated = true;
        arg_bt2_id = -arg_bt2_id;
      }
      if (terms_.find(arg_bt2_id) == terms_.end()) {
        throw PonoException("Missing term for id "
                            + std::to_string(arg_bt2_id));
      }

      Term term_ = terms_.at(arg_bt2_id);
      if (negated) {
        if (term_->get_sort()->get_sort_kind() == BV) {
          term_ = solver_->make_term(BVNot, term_);
        } else {
          term_ = solver_->make_term(Not, term_);
        }
      }
      termargs.push_back(term_);
    }

    /******************************** Handle special cases
     * ********************************/
    if (bt2_line->tag == BTOR2_TAG_state) {
      new_symbol = "state" + to_string(bt2_line->id);
      if (bt2_line->symbol) {
        orig_symbol = bt2_line->symbol;
      } else {
        auto renaming_lookup_pos = state_renaming_table_.find(bt2_line->id);
        if (renaming_lookup_pos != state_renaming_table_.end())
          orig_symbol = renaming_lookup_pos->second;
        else
          orig_symbol = "";
      }
      Term state = ts_.make_statevar(new_symbol, linesort);
      terms_[bt2_line->id] = state;
      statesvec_.push_back(state);
      symbol_map_[new_symbol] = orig_symbol;
      // will be removed from this map if there's a next function for this state
      no_next_states_[num_states] = state;
      id2statenum[bt2_line->id] = num_states;
      num_states++;
    } else if (bt2_line->tag == BTOR2_TAG_input) {
      new_symbol = "input" + to_string(bt2_line->id);
      orig_symbol = bt2_line->symbol ? bt2_line->symbol : "";
      Term input = ts_.make_inputvar(new_symbol, linesort);
      terms_[bt2_line->id] = input;
      inputsvec_.push_back(input);
      symbol_map_[new_symbol] = orig_symbol;
    } else if (bt2_line->tag == BTOR2_TAG_output) {
      new_symbol = "output" + to_string(bt2_line->id);
      orig_symbol = bt2_line->symbol ? bt2_line->symbol : "";
      try {
        ts_.name_term(new_symbol, termargs[0]);
      }
      catch (PonoException & e) {
        new_symbol = "_out_" + new_symbol;
        ts_.name_term(new_symbol, termargs[0]);
      }
      terms_[bt2_line->id] = termargs[0];
      symbol_map_[new_symbol] = orig_symbol;
    } else if (bt2_line->tag == BTOR2_TAG_sort) {
      switch (bt2_line->sort.tag) {
        case BTOR2_TAG_SORT_bitvec: {
          linesort = solver_->make_sort(BV, bt2_line->sort.bitvec.width);
          sorts_[bt2_line->id] = linesort;
          break;
        }
        case BTOR2_TAG_SORT_array: {
          linesort =
              solver_->make_sort(ARRAY,
                                 sorts_.at(bt2_line->sort.array.index),
                                 sorts_.at(bt2_line->sort.array.element));
          sorts_[bt2_line->id] = linesort;
          break;
        }
        default:
          // TODO: maybe only check this in debug? or could always check cause
          // it's really bad
          throw PonoException("Unknown sort tag");
      }
    } else if (bt2_line->tag == BTOR2_TAG_constraint) {
      Term constraint = bv_to_bool(termargs[0]);

      // BTOR2 allows constraints over inputs
      // in Pono these need to be promoted to state variables
      UnorderedTermSet free_vars;
      get_free_symbolic_consts(constraint, free_vars);
      for (const auto & v : free_vars) {
        if (ts_.is_input_var(v)) {
          ts_.promote_inputvar(v);
        }
      }

      ts_.add_constraint(constraint);
      terms_[bt2_line->id] = constraint;
    } else if (bt2_line->tag == BTOR2_TAG_init) {
      if (termargs.size() != 2) {
        throw PonoException("Expecting two term arguments to init");
      } else if (linesort != termargs[0]->get_sort()) {
        throw PonoException(
            "Expecting to init sort to match first argument's sort");
      }

      Term inbt2_iteq;

      if (linesort->get_sort_kind() == BV) {
        inbt2_iteq = solver_->make_term(Equal, termargs);
      } else if (linesort->get_sort_kind() == ARRAY) {
        if (termargs[1]->get_sort()->get_sort_kind() == BV) {
          inbt2_iteq = solver_->make_term(
              Equal, termargs[0], solver_->make_term(termargs[1], linesort));
        } else {
          inbt2_iteq = solver_->make_term(Equal, termargs[0], termargs[1]);
        }
      } else {
        throw PonoException("Unhandled sort: "
                            + termargs[0]->get_sort()->to_string());
      }

      ts_.constrain_init(inbt2_iteq);
      terms_[bt2_line->id] = inbt2_iteq;

    } else if (bt2_line->tag == BTOR2_TAG_next) {
      if (termargs.size() != 2) {
        throw PonoException("Expecting two arguments to next");
      }

      Term t0 = termargs[0];
      Term t1 = termargs[1];
      Sort s0 = t0->get_sort();
      Sort s1 = t1->get_sort();
      SortKind sk0 = s0->get_sort_kind();
      SortKind sk1 = s1->get_sort_kind();

      if (s0 == s1) {
        ts_.assign_next(t0, t1);
        terms_[bt2_line->id] = t1;
      } else if (((sk0 == BV) && (sk1 == BOOL))
                 || ((sk0 == BOOL) && (sk1 == BV))) {
        // need to cast
        ts_.assign_next(bool_to_bv(t0), bool_to_bv(t1));
        terms_[bt2_line->id] = bool_to_bv(t1);
      } else {
        throw PonoException("Got two different sorts in next update.");
      }
      no_next_states_.erase(id2statenum.at(bt2_line->args[0]));
    } else if (bt2_line->tag == BTOR2_TAG_bad) {
      Term bad = bv_to_bool(termargs[0]);
      Term prop = solver_->make_term(Not, bad);

      // BTOR2 allows properties (negation of bad) over inputs
      // in Pono these need to be promoted to state variables
      UnorderedTermSet free_vars;
      get_free_symbolic_consts(bad, free_vars);
      for (const auto & v : free_vars) {
        if (ts_.is_input_var(v)) {
          ts_.promote_inputvar(v);
        }
      }

      propvec_.push_back(prop);
      terms_[bt2_line->id] = prop;
    } else if (bt2_line->tag == BTOR2_TAG_justice) {
      auto & justice = justicevec_.emplace_back();
      justice.reserve(termargs.size());
      std::transform(termargs.begin(),
                     termargs.end(),
                     std::back_inserter(justice),
                     [&](auto t) { return bv_to_bool(t); });
    } else if (bt2_line->tag == BTOR2_TAG_fair) {
      std::cerr << "Warning: ignoring fair term" << std::endl;
      fairvec_.push_back(termargs[0]);
      terms_[bt2_line->id] = termargs[0];
    } else if (bt2_line->constant) {
      terms_[bt2_line->id] = solver_->make_term(
          bt2_line->constant, linesort, basemap.at(bt2_line->tag));
    } else if (bt2_line->tag == BTOR2_TAG_one) {
      terms_[bt2_line->id] = solver_->make_term(1, linesort);
    } else if (bt2_line->tag == BTOR2_TAG_ones) {
      terms_[bt2_line->id] =
          solver_->make_term(string(linesort->get_width(), '1'), linesort, 2);
    } else if (bt2_line->tag == BTOR2_TAG_zero) {
      terms_[bt2_line->id] = solver_->make_term(0, linesort);
    } else if (bt2_line->tag == BTOR2_TAG_eq) {
      if (termargs.size() != 2) {
        throw PonoException("Expecting two arguments to eq");
      }
      Term t0 = termargs[0];
      Term t1 = termargs[1];
      Sort s0 = t0->get_sort();
      Sort s1 = t1->get_sort();
      SortKind sk0 = s0->get_sort_kind();
      SortKind sk1 = s1->get_sort_kind();

      if (s0 != s1) {
        if (((sk0 == BV) && (sk1 == BOOL)) || ((sk0 == BOOL) && (sk1 == BV))) {
          // cast to bit-vectors
          t0 = bool_to_bv(t0);
          t1 = bool_to_bv(t1);
          sk0 = BV;
          sk1 = BV;
        } else {
          throw PonoException(
              "Expecting arguments to eq to have the same sort");
        }
      }

      if (sk0 == BV) {
        terms_[bt2_line->id] = solver_->make_term(BVComp, t0, t1);
      } else {
        terms_[bt2_line->id] = solver_->make_term(Equal, t0, t1);
      }
    } else if (bt2_line->tag == BTOR2_TAG_neq) {
      if (termargs.size() != 2) {
        throw PonoException("Expecting two arguments to neq");
      }
      Term t0 = termargs[0];
      Term t1 = termargs[1];
      Sort s0 = t0->get_sort();
      Sort s1 = t1->get_sort();
      SortKind sk0 = s0->get_sort_kind();
      SortKind sk1 = s1->get_sort_kind();

      if (s0 != s1) {
        if (((sk0 == BV) && (sk1 == BOOL)) || ((sk0 == BOOL) && (sk1 == BV))) {
          // cast to bit-vectors
          t0 = bool_to_bv(t0);
          t1 = bool_to_bv(t1);
          sk0 = BV;
          sk1 = BV;
        } else {
          throw PonoException(
              "Expecting arguments to neq to have the same sort");
        }
      }

      if (sk0 == BV) {
        terms_[bt2_line->id] =
            solver_->make_term(BVNot, solver_->make_term(BVComp, t0, t1));
      } else {
        terms_[bt2_line->id] = solver_->make_term(Distinct, t0, t1);
      }
    } else if (bt2_line->tag == BTOR2_TAG_slice) {
      terms_[bt2_line->id] =
          solver_->make_term(Op(Extract, bt2_line->args[1], bt2_line->args[2]),
                             bool_to_bv(termargs[0]));
    } else if (bt2_line->tag == BTOR2_TAG_sext) {
      terms_[bt2_line->id] = solver_->make_term(
          Op(Sign_Extend, bt2_line->args[1]), bool_to_bv(termargs[0]));
    } else if (bt2_line->tag == BTOR2_TAG_uext) {
      terms_[bt2_line->id] = solver_->make_term(
          Op(Zero_Extend, bt2_line->args[1]), bool_to_bv(termargs[0]));
    } else if (bt2_line->tag == BTOR2_TAG_rol) {
      terms_[bt2_line->id] = solver_->make_term(
          Op(Rotate_Left, bt2_line->args[1]), bool_to_bv(termargs[0]));
    } else if (bt2_line->tag == BTOR2_TAG_ror) {
      terms_[bt2_line->id] = solver_->make_term(
          Op(Rotate_Right, bt2_line->args[1]), bool_to_bv(termargs[0]));
    } else if (bt2_line->tag == BTOR2_TAG_inc) {
      Term t = bool_to_bv(termargs[0]);
      terms_[bt2_line->id] =
          solver_->make_term(BVAdd, t, solver_->make_term(1, t->get_sort()));
    } else if (bt2_line->tag == BTOR2_TAG_dec) {
      Term t = bool_to_bv(termargs[0]);
      terms_[bt2_line->id] =
          solver_->make_term(BVSub, t, solver_->make_term(1, t->get_sort()));
    } else if (bt2_line->tag == BTOR2_TAG_implies) {
      if (termargs.size() != 2) {
        throw PonoException("Expecting two arguments to implies");
      }
      Term t0 = termargs[0];
      Term t1 = termargs[1];
      Sort s0 = t0->get_sort();
      Sort s1 = t1->get_sort();
      SortKind sk0 = s0->get_sort_kind();
      SortKind sk1 = s1->get_sort_kind();

      if (s0 != s1) {
        if (((sk0 == BV) && (sk1 == BOOL)) || ((sk0 == BOOL) && (sk1 == BV))) {
          // cast to bools
          t0 = bv_to_bool(t0);
          t1 = bv_to_bool(t1);
          sk0 = BOOL;
          sk1 = BOOL;
        } else {
          throw PonoException(
              "Expecting arguments to implies to have the same sort");
        }
      }

      if (sk0 == BV) {
        terms_[bt2_line->id] =
            solver_->make_term(BVOr, solver_->make_term(BVNot, t0), t1);
      } else {
        terms_[bt2_line->id] = solver_->make_term(Implies, t0, t1);
      }
    } else if (bt2_line->tag == BTOR2_TAG_redand) {
      Term t = bool_to_bv(termargs[0]);
      Term ones = solver_->make_term(
          std::string(t->get_sort()->get_width(), '1'), t->get_sort(), 2);
      terms_[bt2_line->id] = solver_->make_term(BVComp, t, ones);
    } else if (bt2_line->tag == BTOR2_TAG_redor) {
      Term t = bool_to_bv(termargs[0]);
      Term zero = solver_->make_term(0, t->get_sort());
      terms_[bt2_line->id] = solver_->make_term(Distinct, t, zero);
    } else if (bt2_line->tag == BTOR2_TAG_redxor) {
      Term t = bool_to_bv(termargs[0]);
      unsigned int width = t->get_sort()->get_width();
      Term res = solver_->make_term(Op(Extract, width - 1, width - 1), t);
      for (int i = width - 2; i >= 0; i--) {
        res = solver_->make_term(
            BVXor, res, solver_->make_term(Op(Extract, i, i), t));
      }
      terms_[bt2_line->id] = res;
    } else if (bt2_line->tag == BTOR2_TAG_ite) {
      Term cond = bv_to_bool(termargs[0]);
      // Always cast to bit-vectors because mathsat doesn't support ite over
      // bools
      Term t1 = bool_to_bv(termargs[1]);
      Term t2 = bool_to_bv(termargs[2]);
      terms_[bt2_line->id] = solver_->make_term(Ite, cond, t1, t2);
    } else if (bt2_line->tag == BTOR2_TAG_uaddo) {
      Term t0 = bool_to_bv(termargs[0]);
      Term t1 = bool_to_bv(termargs[1]);

      int orig_width = t0->get_sort()->get_width();

      t0 = solver_->make_term(Op(Zero_Extend, 1), t0);
      t1 = solver_->make_term(Op(Zero_Extend, 1), t1);

      Term sum = solver_->make_term(BVAdd, t0, t1);
      // overflow occurs if there's a carry out bit
      terms_[bt2_line->id] =
          solver_->make_term(Op(Extract, orig_width, orig_width), sum);
    } else if (bt2_line->tag == BTOR2_TAG_saddo) {
      // From https://www.doc.ic.ac.uk/~eedwards/compsys/arithmetic/index.html
      Term t0 = bool_to_bv(termargs[0]);
      Term t1 = bool_to_bv(termargs[1]);
      int width = t0->get_sort()->get_width();
      Term sum = solver_->make_term(BVAdd, t0, t1);
      // overflow occurs if
      // both operands are positive and the result is negative or
      // both operands are negative and the result is positive
      Term t0_top = solver_->make_term(Op(Extract, width - 1, width - 1), t0);
      Term t1_top = solver_->make_term(Op(Extract, width - 1, width - 1), t1);
      Term sum_top = solver_->make_term(Op(Extract, width - 1, width - 1), sum);
      terms_[bt2_line->id] =
          solver_->make_term(Equal,
                             solver_->make_term(Equal, t0_top, t1_top),
                             solver_->make_term(Distinct, t0_top, sum_top));
    } else if (bt2_line->tag == BTOR2_TAG_sdivo) {
      Term t0 = bool_to_bv(termargs[0]);
      Term t1 = bool_to_bv(termargs[1]);
      int width = t0->get_sort()->get_width();
      Sort sort = solver_->make_sort(BV, width);
      Term sum = solver_->make_term(BVAdd, t0, t1);
      // overflow occurs if
      // t0 is int_min (e.g. 1 followed by zeros) and
      // t1 is -1
      std::string sint_min("1");
      sint_min += std::string(width - 1, '0');
      std::string snegone = std::string(width, '1');
      Term int_min = solver_->make_term(sint_min, sort, 2);
      Term negone = solver_->make_term(snegone, sort, 2);
      terms_[bt2_line->id] =
          solver_->make_term(And,
                             solver_->make_term(Equal, t0, int_min),
                             solver_->make_term(Equal, t1, negone));
    } else if (bt2_line->tag == BTOR2_TAG_umulo) {
      // from Hacker's Delight
      // overflow if hi(x*y) != 0
      Term t0 = bool_to_bv(termargs[0]);
      Term t1 = bool_to_bv(termargs[1]);

      int orig_width = t0->get_sort()->get_width();

      t0 = solver_->make_term(Op(Zero_Extend, orig_width), t0);
      t1 = solver_->make_term(Op(Zero_Extend, orig_width), t1);

      Term prod = solver_->make_term(BVMul, t0, t1);
      // overflow occurs if the upper bits are non-zero
      terms_[bt2_line->id] = solver_->make_term(
          Distinct,
          solver_->make_term(Op(Extract, 2 * orig_width - 1, orig_width), prod),
          solver_->make_term(0, solver_->make_sort(BV, orig_width)));
    } else if (bt2_line->tag == BTOR2_TAG_smulo) {
      // from Hacker's Delight
      // overflow if hi(x*y) != (lo(x*y) >>s (width-1))
      Term t0 = bool_to_bv(termargs[0]);
      Term t1 = bool_to_bv(termargs[1]);

      int orig_width = t0->get_sort()->get_width();

      t0 = solver_->make_term(Op(Zero_Extend, orig_width), t0);
      t1 = solver_->make_term(Op(Zero_Extend, orig_width), t1);

      Term prod = solver_->make_term(BVMul, t0, t1);
      Term hi =
          solver_->make_term(Op(Extract, 2 * orig_width - 1, orig_width), prod);
      Term lo = solver_->make_term(Op(Extract, orig_width - 1, 0), prod);
      terms_[bt2_line->id] = solver_->make_term(
          Distinct,
          hi,
          solver_->make_term(
              BVAshr,
              lo,
              solver_->make_term(orig_width - 1,
                                 solver_->make_sort(BV, orig_width))));
    } else if (bt2_line->tag == BTOR2_TAG_usubo) {
      // From
      // https://github.com/Boolector/boolector/blob/cd757d099433d95ffdb2a839504b220eff18ee51/src/btorexp.c#L1236
      Term t0 = bool_to_bv(termargs[0]);
      Term t1 = bool_to_bv(termargs[1]);
      unsigned int width = t0->get_sort()->get_width();
      Sort sort = solver_->make_sort(BV, width + 1);
      t0 = solver_->make_term(Op(Zero_Extend, 1), t0);
      t1 = solver_->make_term(Op(Zero_Extend, 1), t1);
      Term one = solver_->make_term(1, sort);
      Term add1 = solver_->make_term(BVAdd, t1, one);
      Term add2 = solver_->make_term(BVAdd, t0, add1);
      terms_[bt2_line->id] = solver_->make_term(
          BVNot, solver_->make_term(Op(Extract, width, width), add2));
    } else if (bt2_line->tag == BTOR2_TAG_ssubo) {
      // From https://www.doc.ic.ac.uk/~eedwards/compsys/arithmetic/index.html
      // overflow occurs if signs are different and subtrahend sign matches
      // result sign
      // TODO: check this, not sure if it can overflow when signs aren't
      // different
      //       seems like maybe not but should double check
      Term t0 = bool_to_bv(termargs[0]);
      Term t1 = bool_to_bv(termargs[1]);

      int width = t0->get_sort()->get_width();

      Term diff = solver_->make_term(BVSub, t0, t1);
      Term t0_top = solver_->make_term(Op(Extract, width - 1, width - 1), t0);
      Term t1_top = solver_->make_term(Op(Extract, width - 1, width - 1), t1);
      Term diff_top =
          solver_->make_term(Op(Extract, width - 1, width - 1), diff);
      terms_[bt2_line->id] =
          solver_->make_term(And,
                             solver_->make_term(Distinct, t0_top, t1_top),
                             solver_->make_term(Equal, t1_top, diff_top));
    } else if (bt2_line->tag == BTOR2_TAG_read) {
      Term arr = termargs[0];
      Term idx = bool_to_bv(termargs[1]);
      terms_[bt2_line->id] = solver_->make_term(Select, arr, idx);
    } else if (bt2_line->tag == BTOR2_TAG_write) {
      Term arr = termargs[0];
      Term idx = bool_to_bv(termargs[1]);
      Term elem = bool_to_bv(termargs[2]);
      terms_[bt2_line->id] = solver_->make_term(Store, arr, idx, elem);
    }
    /******************************** Handle general case
     ********************************/
    else {
      if (!termargs.size()) {
        throw PonoException("Expecting non-zero number of terms");
      }

      if (boolopmap.find(bt2_line->tag) != boolopmap.end()) {
        // TODO: potentially remove this, have to treat specially for boolector
        // vs other solvers anyway
        if (bvopmap.find(bt2_line->tag) == bvopmap.end()) {
          // only a boolean op
          // convert all to bools
          for (size_t i = 0; i < termargs.size(); i++) {
            termargs[i] = bv_to_bool(termargs[i]);
          }
        } else {
          termargs = lazy_convert(termargs);
        }

        SortKind sk = termargs[0]->get_sort()->get_sort_kind();
        if (sk == BV) {
          terms_[bt2_line->id] =
              solver_->make_term(bvopmap.at(bt2_line->tag), termargs);
        } else if (sk == BOOL) {
          terms_[bt2_line->id] =
              solver_->make_term(boolopmap.at(bt2_line->tag), termargs);
        } else {
          throw PonoException("Unexpected sort");
        }
      } else {
        for (int i = 0; i < termargs.size(); i++) {
          termargs[i] = bool_to_bv(termargs[i]);
        }
        terms_[bt2_line->id] =
            solver_->make_term(bvopmap.at(bt2_line->tag), termargs);
      }
    }

    // use the symbol to name the term (if applicable)
    // input, output, and state already named
    if (bt2_line->symbol && bt2_line->tag != BTOR2_TAG_input
        && bt2_line->tag != BTOR2_TAG_output && bt2_line->tag != BTOR2_TAG_state
        && terms_.find(bt2_line->id) != terms_.end()) {
      try {
        new_symbol = "term" + to_string(bt2_line->id);
        orig_symbol = bt2_line->symbol;
        ts_.name_term(new_symbol, terms_.at(bt2_line->id));
        symbol_map_[new_symbol] = orig_symbol;
      }
      catch (PonoException & e) {
        logger.log(1, "BTOR2Encoder Warning: {}", e.what());
      }
    }

    // sort and justice tags should be the only that don't populate terms_
    assert(bt2_line->tag == BTOR2_TAG_sort || bt2_line->tag == BTOR2_TAG_justice
           || terms_.find(bt2_line->id) != terms_.end());
  }

  fclose(input_file);
  btor2parser_delete(reader);
}
}  // namespace pono
