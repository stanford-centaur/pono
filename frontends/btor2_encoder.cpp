#include "btor2_encoder.h"

#include <iostream>

using namespace smt;
using namespace std;

namespace cosa {

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
    { BTOR2_TAG_neq, Distinct },
    { BTOR2_TAG_neg, BVNeg },
    //{ BTOR2_TAG_next, },
    { BTOR2_TAG_nor, BVNor },
    { BTOR2_TAG_not, BVNot },
    //{ BTOR2_TAG_one, },
    //{ BTOR2_TAG_ones, },
    { BTOR2_TAG_or, BVOr },
    //{ BTOR2_TAG_output, },
    // { BTOR2_TAG_read, Select }, // handle specially -- make sure it's casted to bv
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
    //{ BTOR2_TAG_write, Store }, // handle specially -- make sure it's casted to bv
    { BTOR2_TAG_xnor, BVXnor },
    { BTOR2_TAG_xor, BVXor },
    //{ BTOR2_TAG_zero, }
});

const unordered_map<Btor2Tag, smt::PrimOp> boolopmap(
    { { BTOR2_TAG_and, And },
      { BTOR2_TAG_or, Or },
      { BTOR2_TAG_xor, Xor },
      { BTOR2_TAG_not, Not },
      // { BTOR2_TAG_implies, Implies },
      { BTOR2_TAG_iff, Iff },
      // { BTOR2_TAG_eq, Equal }, // handling specially -- could have array
      // arguments
      { BTOR2_TAG_neq, Distinct } });

Term BTOR2Encoder::bool_to_bv(const Term & t) const
{
  if (t->get_sort()->get_sort_kind() == BOOL)
  {
    Sort bv1sort = solver_->make_sort(BV, 1);
    return solver_->make_term(
        Ite, t, solver_->make_term(1, bv1sort), solver_->make_term(0, bv1sort));
  }
  else
  {
    return t;
  }
}

Term BTOR2Encoder::bv_to_bool(const Term & t) const
{
  Sort sort = t->get_sort();
  if (sort->get_sort_kind() == BV)
  {
    if (sort->get_width() != 1)
    {
      throw CosaException("Can't convert non-width 1 bitvector to bool.");
    }
    return solver_->make_term(
        Equal, t, solver_->make_term(1, solver_->make_sort(BV, 1)));
  }
  else
  {
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
  for (auto t : tvec)
  {
    res.push_back(t);

    sort = t->get_sort();
    sortset.insert(sort);

    if (sort->get_sort_kind() == BOOL)
    {
      num_bools++;
    }
    else if (!(sort->get_sort_kind() == BV && sort->get_width() == 1))
    {
      wide_bvs = true;
    }
  }

  if (sortset.size() > 1)
  {
    if (num_bools > tvec.size() / 2 && !wide_bvs)
    {
      for (size_t i = 0; i < res.size(); i++)
      {
        res[i] = bv_to_bool(res[i]);
      }
    }
    else
    {
      for (size_t i = 0; i < res.size(); i++)
      {
        res[i] = bool_to_bv(res[i]);
      }
    }
  }

  return res;
}

void BTOR2Encoder::parse(const std::string filename)
{
  FILE * input_file = fopen(filename.c_str(), "r");

  if (!input_file)
  {
    throw CosaException("Could not open " + filename);
  }

  reader_ = btor2parser_new();

  if (!btor2parser_read_lines(reader_, input_file))
  {
    throw CosaException("Error parsing btor file.");
  }

  uint64_t num_states = 0;
  std::unordered_map<int64_t, uint64_t> id2statenum;

  it_ = btor2parser_iter_init(reader_);
  while ((l_ = btor2parser_iter_next(&it_)))
  {
    /******************************** Identify sort
     * ********************************/
    if (l_->tag != BTOR2_TAG_sort && l_->sort.id)
    {
      linesort_ = sorts_.at(l_->sort.id);
    }

    /******************************** Gather term arguments
     * ********************************/
    termargs_.clear();
    termargs_.reserve(l_->nargs);
    for (i_ = 0; i_ < l_->nargs; i_++)
    {
      negated_ = false;
      idx_ = l_->args[i_];
      if (idx_ < 0)
      {
        negated_ = true;
        idx_ = -idx_;
      }
      if (terms_.find(idx_) == terms_.end())
      {
        throw CosaException("Missing term for id " + std::to_string(idx_));
      }

      Term term_ = terms_.at(idx_);
      if (negated_)
      {
        if (term_->get_sort()->get_sort_kind() == BV)
        {
          term_ = solver_->make_term(BVNot, term_);
        }
        else
        {
          term_ = solver_->make_term(Not, term_);
        }
      }
      termargs_.push_back(term_);
    }

    /******************************** Handle special cases
     * ********************************/
    if (l_->tag == BTOR2_TAG_state)
    {
      if (l_->symbol)
      {
        symbol_ = l_->symbol;
      }
      else
      {
        symbol_ = "state" + to_string(l_->id);
      }

      Term state = rts_.make_state(symbol_, linesort_);
      terms_[l_->id] = state;
      statesvec_.push_back(state);
      // will be removed from this map if there's a next function for this state
      no_next_states_[num_states] = state;
      id2statenum[l_->id] = num_states;
      num_states++;
    }
    else if (l_->tag == BTOR2_TAG_input)
    {
      if (l_->symbol)
      {
        symbol_ = l_->symbol;
      }
      else
      {
        symbol_ = "input" + to_string(l_->id);
      }
      Term input = rts_.make_input(symbol_, linesort_);
      terms_[l_->id] = input;
      inputsvec_.push_back(input);
    }
    else if (l_->tag == BTOR2_TAG_output)
    {
      if (l_->symbol)
      {
        symbol_ = l_->symbol;
      }
      else
      {
        symbol_ = "output" + to_string(l_->id);
      }

      rts_.name_term(symbol_, termargs_[0]);
    }
    else if (l_->tag == BTOR2_TAG_sort)
    {
      switch (l_->sort.tag)
      {
        case BTOR2_TAG_SORT_bitvec:
        {
          linesort_ = solver_->make_sort(BV, l_->sort.bitvec.width);
          sorts_[l_->id] = linesort_;
          break;
        }
        case BTOR2_TAG_SORT_array:
        {
          linesort_ = solver_->make_sort(ARRAY,
                                         sorts_.at(l_->sort.array.index),
                                         sorts_.at(l_->sort.array.element));
          sorts_[l_->id] = linesort_;
          break;
        }
        default:
          // TODO: maybe only check this in debug? or could always check cause
          // it's really bad
          throw CosaException("Unknown sort tag");
      }
    }
    else if (l_->tag == BTOR2_TAG_constraint)
    {
      rts_.add_invar(bv_to_bool(termargs_[0]));
    }
    else if (l_->tag == BTOR2_TAG_init)
    {
      if (termargs_.size() != 2)
      {
        throw CosaException("Expecting two term arguments to init");
      }
      else if (linesort_ != termargs_[0]->get_sort())
      {
        throw CosaException(
            "Expecting to init sort to match first argument's sort");
      }

      if (linesort_->get_sort_kind() == BV)
      {
        rts_.constrain_init(solver_->make_term(Equal, termargs_));
      }
      else if (linesort_->get_sort_kind() == ARRAY)
      {
        if (termargs_[1]->get_sort()->get_sort_kind() == BV)
        {
          rts_.constrain_init(
              solver_->make_term(Equal,
                                 termargs_[0],
                                 solver_->make_term(termargs_[1], linesort_)));
        }
        else
        {
          rts_.constrain_init(
              solver_->make_term(Equal, termargs_[0], termargs_[1]));
        }
      }
      else
      {
        throw CosaException("Unhandled sort: "
                            + termargs_[0]->get_sort()->to_string());
      }
    }
    else if (l_->tag == BTOR2_TAG_next)
    {
      if (termargs_.size() != 2)
      {
        throw CosaException("Expecting two arguments to next");
      }

      Term t0 = termargs_[0];
      Term t1 = termargs_[1];
      Sort s0 = t0->get_sort();
      Sort s1 = t1->get_sort();
      SortKind sk0 = s0->get_sort_kind();
      SortKind sk1 = s1->get_sort_kind();

      if (s0 == s1)
      {
        rts_.set_next(t0, t1);
      }
      else if (((sk0 == BV) && (sk1 == BOOL)) || ((sk0 == BOOL) && (sk1 == BV)))
      {
        // need to cast
        rts_.set_next(bool_to_bv(t0), bool_to_bv(t1));
      }
      else
      {
        throw CosaException("Got two different sorts in next update.");
      }
      no_next_states_.erase(id2statenum.at(l_->args[0]));
    }
    else if (l_->tag == BTOR2_TAG_bad)
    {
      badvec_.push_back(bv_to_bool(termargs_[0]));
    }
    else if (l_->tag == BTOR2_TAG_justice)
    {
      std::cout << "Warning: ignoring justice term" << std::endl;
      justicevec_.push_back(termargs_[0]);
    }
    else if (l_->tag == BTOR2_TAG_fair)
    {
      std::cout << "Warning: ignoring fair term" << std::endl;
      fairvec_.push_back(termargs_[0]);
    }
    else if (l_->constant)
    {
      terms_[l_->id] =
          solver_->make_term(l_->constant, linesort_, basemap.at(l_->tag));
    }
    else if (l_->tag == BTOR2_TAG_one)
    {
      terms_[l_->id] = solver_->make_term(1, linesort_);
    }
    else if (l_->tag == BTOR2_TAG_ones)
    {
      terms_[l_->id] =
          solver_->make_term(string(linesort_->get_width(), '1'), linesort_, 2);
    }
    else if (l_->tag == BTOR2_TAG_zero)
    {
      terms_[l_->id] = solver_->make_term(0, linesort_);
    }
    else if (l_->tag == BTOR2_TAG_slice)
    {
      terms_[l_->id] = solver_->make_term(Op(Extract, l_->args[1], l_->args[2]),
                                          bool_to_bv(termargs_[0]));
    }
    else if (l_->tag == BTOR2_TAG_sext)
    {
      terms_[l_->id] = solver_->make_term(Op(Sign_Extend, l_->args[1]),
                                          bool_to_bv(termargs_[0]));
    }
    else if (l_->tag == BTOR2_TAG_uext)
    {
      terms_[l_->id] = solver_->make_term(Op(Zero_Extend, l_->args[1]),
                                          bool_to_bv(termargs_[0]));
    }
    else if (l_->tag == BTOR2_TAG_rol)
    {
      terms_[l_->id] = solver_->make_term(Op(Rotate_Left, l_->args[1]),
                                          bool_to_bv(termargs_[0]));
    }
    else if (l_->tag == BTOR2_TAG_ror)
    {
      terms_[l_->id] = solver_->make_term(Op(Rotate_Right, l_->args[1]),
                                          bool_to_bv(termargs_[0]));
    }
    else if (l_->tag == BTOR2_TAG_inc)
    {
      Term t = bool_to_bv(termargs_[0]);
      terms_[l_->id] =
          solver_->make_term(BVAdd, t, solver_->make_term(1, t->get_sort()));
    }
    else if (l_->tag == BTOR2_TAG_dec)
    {
      Term t = bool_to_bv(termargs_[0]);
      terms_[l_->id] =
          solver_->make_term(BVSub, t, solver_->make_term(1, t->get_sort()));
    }
    else if (l_->tag == BTOR2_TAG_eq)
    {
      if (termargs_.size() != 2)
      {
        throw CosaException("Expecting two arguments to eq");
      }
      Term t0 = termargs_[0];
      Term t1 = termargs_[1];
      Sort s0 = t0->get_sort();
      Sort s1 = t1->get_sort();
      SortKind sk0 = s0->get_sort_kind();
      SortKind sk1 = s1->get_sort_kind();

      if (s0 != s1)
      {
        if (((sk0 == BV) && (sk1 == BOOL)) || ((sk0 == BOOL) && (sk1 == BV)))
        {
          // cast to bit-vectors
          t0 = bool_to_bv(t0);
          t1 = bool_to_bv(t1);
          sk0 = BV;
          sk1 = BV;
        }
        else
        {
          throw CosaException(
              "Expecting arguments to eq to have the same sort");
        }
      }

      if (sk0 == BV)
      {
        terms_[l_->id] = solver_->make_term(BVComp, t0, t1);
      }
      else
      {
        terms_[l_->id] = solver_->make_term(Equal, t0, t1);
      }
    }
    else if (l_->tag == BTOR2_TAG_implies)
    {
      if (termargs_.size() != 2)
      {
        throw CosaException("Expecting two arguments to implies");
      }
      Term t0 = termargs_[0];
      Term t1 = termargs_[1];
      Sort s0 = t0->get_sort();
      Sort s1 = t1->get_sort();
      SortKind sk0 = s0->get_sort_kind();
      SortKind sk1 = s1->get_sort_kind();

      if (s0 != s1)
      {
        if (((sk0 == BV) && (sk1 == BOOL)) || ((sk0 == BOOL) && (sk1 == BV)))
        {
          // cast to bools
          t0 = bv_to_bool(t0);
          t1 = bv_to_bool(t1);
          sk0 = BOOL;
          sk1 = BOOL;
        }
        else
        {
          throw CosaException(
              "Expecting arguments to implies to have the same sort");
        }
      }

      if (sk0 == BV)
      {
        terms_[l_->id] =
            solver_->make_term(BVOr, solver_->make_term(BVNot, t0), t1);
      }
      else
      {
        terms_[l_->id] = solver_->make_term(Implies, t0, t1);
      }
    }
    else if (l_->tag == BTOR2_TAG_redand)
    {
      Term t = bool_to_bv(termargs_[0]);
      Term ones = solver_->make_term(
          std::string(t->get_sort()->get_width(), '1'), t->get_sort(), 2);
      terms_[l_->id] = solver_->make_term(BVComp, t, ones);
    }
    else if (l_->tag == BTOR2_TAG_redor)
    {
      Term t = bool_to_bv(termargs_[0]);
      Term zero = solver_->make_term(0, t->get_sort());
      terms_[l_->id] = solver_->make_term(Distinct, t, zero);
    }
    else if (l_->tag == BTOR2_TAG_redxor)
    {
      Term t = bool_to_bv(termargs_[0]);
      unsigned int width = t->get_sort()->get_width();
      Term res = solver_->make_term(Op(Extract, width - 1, width - 1), t);
      for (int i = width - 2; i >= 0; i--)
      {
        res = solver_->make_term(
            BVXor, res, solver_->make_term(Op(Extract, i, i), t));
      }
      terms_[l_->id] = res;
    }
    else if (l_->tag == BTOR2_TAG_ite)
    {
      Term cond = bv_to_bool(termargs_[0]);
      // Always cast to bit-vectors because mathsat doesn't support ite over bools
      Term t1 = bool_to_bv(termargs_[1]);
      Term t2 = bool_to_bv(termargs_[2]);
      terms_[l_->id] = solver_->make_term(Ite, cond, t1, t2);
    }
    else if (l_->tag == BTOR2_TAG_uaddo)
    {
      Term t0 = bool_to_bv(termargs_[0]);
      Term t1 = bool_to_bv(termargs_[1]);

      int orig_width = t0->get_sort()->get_width();

      t0 = solver_->make_term(Op(Zero_Extend, 1), termargs_[0]);
      t1 = solver_->make_term(Op(Zero_Extend, 1), termargs_[0]);

      Term sum = solver_->make_term(BVAdd, t0, t1);
      // overflow occurs if there's a carry out bit
      terms_[l_->id] =
          solver_->make_term(Op(Extract, orig_width, orig_width), sum);
    }
    else if (l_->tag == BTOR2_TAG_saddo)
    {
      // From https://www.doc.ic.ac.uk/~eedwards/compsys/arithmetic/index.html
      Term t0 = bool_to_bv(termargs_[0]);
      Term t1 = bool_to_bv(termargs_[1]);
      int width = t0->get_sort()->get_width();
      Term sum = solver_->make_term(BVAdd, t0, t1);
      // overflow occurs if
      // both operands are positive and the result is negative or
      // both operands are negative and the result is positive
      Term t0_top = solver_->make_term(Op(Extract, width - 1, width - 1), t0);
      Term t1_top = solver_->make_term(Op(Extract, width - 1, width - 1), t1);
      Term sum_top = solver_->make_term(Op(Extract, width - 1, width - 1), sum);
      terms_[l_->id] =
          solver_->make_term(Equal,
                             solver_->make_term(Equal, t0_top, t1_top),
                             solver_->make_term(Distinct, t0_top, sum_top));
    }
    else if (l_->tag == BTOR2_TAG_sdivo)
    {
      Term t0 = bool_to_bv(termargs_[0]);
      Term t1 = bool_to_bv(termargs_[1]);
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
      terms_[l_->id] =
          solver_->make_term(And,
                             solver_->make_term(Equal, t0, int_min),
                             solver_->make_term(Equal, t1, negone));
    }
    else if (l_->tag == BTOR2_TAG_umulo)
    {
      // from Hacker's Delight
      // overflow if hi(x*y) != 0
      Term t0 = bool_to_bv(termargs_[0]);
      Term t1 = bool_to_bv(termargs_[1]);

      int orig_width = t0->get_sort()->get_width();

      t0 = solver_->make_term(Op(Zero_Extend, orig_width), termargs_[0]);
      t1 = solver_->make_term(Op(Zero_Extend, orig_width), termargs_[0]);

      Term prod = solver_->make_term(BVMul, t0, t1);
      // overflow occurs if the upper bits are non-zero
      terms_[l_->id] = solver_->make_term(
          Distinct,
          solver_->make_term(Op(Extract, 2 * orig_width - 1, orig_width), prod),
          solver_->make_term(0, solver_->make_sort(BV, orig_width)));
    }
    else if (l_->tag == BTOR2_TAG_smulo)
    {
      // from Hacker's Delight
      // overflow if hi(x*y) != (lo(x*y) >>s (width-1))
      Term t0 = bool_to_bv(termargs_[0]);
      Term t1 = bool_to_bv(termargs_[1]);

      int orig_width = t0->get_sort()->get_width();

      t0 = solver_->make_term(Op(Zero_Extend, orig_width), termargs_[0]);
      t1 = solver_->make_term(Op(Zero_Extend, orig_width), termargs_[0]);

      Term prod = solver_->make_term(BVMul, t0, t1);
      Term hi =
          solver_->make_term(Op(Extract, 2 * orig_width - 1, orig_width), prod);
      Term lo = solver_->make_term(Op(Extract, orig_width - 1, 0), prod);
      terms_[l_->id] = solver_->make_term(
          Distinct,
          hi,
          solver_->make_term(
              BVAshr,
              lo,
              solver_->make_term(orig_width - 1,
                                 solver_->make_sort(BV, orig_width))));
    }
    else if (l_->tag == BTOR2_TAG_usubo)
    {
      // From
      // https://github.com/Boolector/boolector/blob/cd757d099433d95ffdb2a839504b220eff18ee51/src/btorexp.c#L1236
      Term t0 = bool_to_bv(termargs_[0]);
      Term t1 = bool_to_bv(termargs_[1]);
      unsigned int width = t0->get_sort()->get_width();
      Sort sort = solver_->make_sort(BV, width + 1);
      t0 = solver_->make_term(Op(Zero_Extend, 1), t0);
      t1 = solver_->make_term(Op(Zero_Extend, 1), t1);
      Term one = solver_->make_term(1, sort);
      Term add1 = solver_->make_term(BVAdd, t1, one);
      Term add2 = solver_->make_term(BVAdd, t0, add1);
      terms_[l_->id] = solver_->make_term(
          BVNot, solver_->make_term(Op(Extract, width, width), add2));
    }
    else if (l_->tag == BTOR2_TAG_ssubo)
    {
      // From https://www.doc.ic.ac.uk/~eedwards/compsys/arithmetic/index.html
      // overflow occurs if signs are different and subtrahend sign matches
      // result sign
      // TODO: check this, not sure if it can overflow when signs aren't
      // different
      //       seems like maybe not but should double check
      Term t0 = bool_to_bv(termargs_[0]);
      Term t1 = bool_to_bv(termargs_[1]);

      int width = t0->get_sort()->get_width();

      Term diff = solver_->make_term(BVSub, t0, t1);
      Term t0_top = solver_->make_term(Op(Extract, width - 1, width - 1), t0);
      Term t1_top = solver_->make_term(Op(Extract, width - 1, width - 1), t1);
      Term diff_top =
          solver_->make_term(Op(Extract, width - 1, width - 1), diff);
      terms_[l_->id] =
          solver_->make_term(And,
                             solver_->make_term(Distinct, t0_top, t1_top),
                             solver_->make_term(Equal, t1_top, diff_top));
    }
    else if (l_->tag == BTOR2_TAG_read)
    {
      Term arr = termargs_[0];
      Term idx = bool_to_bv(termargs_[1]);
      terms_[l_->id] = solver_->make_term(Select, arr, idx);
    }
    else if (l_->tag == BTOR2_TAG_write)
    {
      Term arr = termargs_[0];
      Term idx = bool_to_bv(termargs_[1]);
      Term elem = bool_to_bv(termargs_[2]);
      terms_[l_->id] = solver_->make_term(Store, arr, idx, elem);
    }
    /******************************** Handle general case
     ********************************/
    else
    {
      if (!termargs_.size())
      {
        throw CosaException("Expecting non-zero number of terms");
      }

      if (boolopmap.find(l_->tag) != boolopmap.end())
      {
        // TODO: potentially remove this, have to treat specially for boolector
        // vs other solvers anyway
        if (bvopmap.find(l_->tag) == bvopmap.end())
        {
          // only a boolean op
          // convert all to bools
          for (size_t i = 0; i < termargs_.size(); i++)
          {
            termargs_[i] = bv_to_bool(termargs_[i]);
          }
        }
        else
        {
          termargs_ = lazy_convert(termargs_);
        }

        SortKind sk = termargs_[0]->get_sort()->get_sort_kind();
        if (sk == BV)
        {
          terms_[l_->id] = solver_->make_term(bvopmap.at(l_->tag), termargs_);
        }
        else if (sk == BOOL)
        {
          terms_[l_->id] = solver_->make_term(boolopmap.at(l_->tag), termargs_);
        }
        else
        {
          throw CosaException("Unexpected sort");
        }
      }
      else
      {
        for (int i = 0; i < termargs_.size(); i++)
        {
          termargs_[i] = bool_to_bv(termargs_[i]);
        }
        terms_[l_->id] = solver_->make_term(bvopmap.at(l_->tag), termargs_);
      }
    }
  }

  fclose(input_file);
  btor2parser_delete(reader_);
}
}  // namespace cosa
