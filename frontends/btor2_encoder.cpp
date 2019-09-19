#include "btor2_encoder.h"

#include <iostream>

using namespace smt;
using namespace std;

namespace cosa
{

// Maps for use in conversion
const unordered_map<Btor2Tag, int> basemap(
                                           { {BTOR2_TAG_const, 2},
                                             {BTOR2_TAG_constd, 10},
                                             {BTOR2_TAG_consth, 16} });
const unordered_set<Btor2Tag> overflow_ops({BTOR2_TAG_uaddo,
                                            BTOR2_TAG_saddo,
                                            //BTOR2_TAG_udivo, // not supposed to exist?
                                            BTOR2_TAG_sdivo,
                                            BTOR2_TAG_umulo,
                                            BTOR2_TAG_smulo,
                                            BTOR2_TAG_usubo,
                                            BTOR2_TAG_ssubo});

const unordered_map<Btor2Tag, smt::PrimOp> bvopmap(
                                              {  { BTOR2_TAG_add, BVAdd },
                                                 { BTOR2_TAG_and, BVAnd },
                                                 // { BTOR2_TAG_bad, },
                                                 { BTOR2_TAG_concat, Concat },
                                                 //{ BTOR2_TAG_const, },
                                                 //{ BTOR2_TAG_constraint, },
                                                 //{ BTOR2_TAG_constd, },
                                                 //{ BTOR2_TAG_consth, },
                                                 //{ BTOR2_TAG_dec, },
                                                 { BTOR2_TAG_eq, BVComp },
                                                 //{ BTOR2_TAG_fair, },
                                                 { BTOR2_TAG_iff, Iff },
                                                 { BTOR2_TAG_implies, Implies },
                                                 //{ BTOR2_TAG_inc, },
                                                 //{ BTOR2_TAG_init, },
                                                 //{ BTOR2_TAG_input, },
                                                 //{ BTOR2_TAG_ite, Ite },
                                                 //{ BTOR2_TAG_justice, },
                                                 { BTOR2_TAG_mul, BVMul },
                                                 { BTOR2_TAG_nand, BVNand },
                                                 { BTOR2_TAG_neq, Distinct},
                                                 { BTOR2_TAG_neg, BVNeg },
                                                 //{ BTOR2_TAG_next, },
                                                 { BTOR2_TAG_nor, BVNor },
                                                 { BTOR2_TAG_not, BVNot },
                                                 //{ BTOR2_TAG_one, },
                                                 //{ BTOR2_TAG_ones, },
                                                 { BTOR2_TAG_or, BVOr },
                                                 //{ BTOR2_TAG_output, },
                                                 { BTOR2_TAG_read, Select },
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
                                                 { BTOR2_TAG_write, Store },
                                                 { BTOR2_TAG_xnor, BVXnor },
                                                 { BTOR2_TAG_xor, BVXor },
                                                 //{ BTOR2_TAG_zero, }
                                              });

const unordered_map<Btor2Tag, smt::PrimOp> boolopmap(
                                                  { { BTOR2_TAG_and, And },
                                                    { BTOR2_TAG_or, Or },
                                                    { BTOR2_TAG_xor, Xor },
                                                    { BTOR2_TAG_not, Not },
                                                    { BTOR2_TAG_implies, Implies },
                                                    { BTOR2_TAG_iff, Iff },
                                                    { BTOR2_TAG_eq, Equal },
                                                    { BTOR2_TAG_neq, Distinct }
                                                  });

Term BTOR2Encoder::bool_to_bv(Term t)
{
  if (t->get_sort()->get_sort_kind() == BOOL)
  {
    Sort bv1sort = solver_->make_sort(BV, 1);
    return solver_->make_term(Ite, solver_->make_value(1, bv1sort), solver_->make_value(0, bv1sort));
  }
  else
  {
    return t;
  }
}

Term BTOR2Encoder::bv_to_bool(Term t)
{
  Sort sort = t->get_sort();
  if (sort->get_sort_kind() == BV)
  {
    if (sort->get_width() != 1)
    {
      throw CosaException("Can't convert non-width 1 bitvector to bool.");
    }
    return solver_->make_term(Equal, t, solver_->make_value(1, solver_->make_sort(BV, 1)));
  }
  else
  {
    return t;
  }
}

TermVec BTOR2Encoder::lazy_convert(TermVec tvec)
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
    else if (!(sort->get_sort_kind() == BV &&
               sort->get_width() == 1))
    {
      wide_bvs = true;
    }
  }

  if (sortset.size() > 1)
  {
    if (num_bools > tvec.size()/2 && !wide_bvs)
    {
      for(auto t : tvec)
      {
        res.push_back(bv_to_bool(t));
      }
    }
    else
    {
      for(auto t : tvec)
      {
        res.push_back(bool_to_bv(t));
      }
    }
  }

  return res;
}

void BTOR2Encoder::parse(std::string filename)
{
  FILE* input_file = fopen(filename.c_str(), "r");

  if(!input_file)
  {
    throw CosaException("Could not open " + filename);
  }

  reader_ = btor2parser_new();

  if(!btor2parser_read_lines(reader_, input_file))
  {
    throw CosaException("Error parsing btor file.");
  }

  it_ = btor2parser_iter_init(reader_);
  while((l_ = btor2parser_iter_next(&it_)))
  {
    /******************************** Identify sort ********************************/
    if (l_->tag != BTOR2_TAG_sort && l_->sort.id)
    {
      linesort_=sorts_.at(l_->sort.id);
    }

    /******************************** Gather term arguments ********************************/
    termargs_.clear();
    termargs_.reserve(l_->nargs);
    for(i_ = 0; i_ < l_->nargs; i_++)
    {
      negated_ = false;
      idx_ = l_->args[i_];
      if (idx_ < 0) {
        negated_ = true;
        idx_ = -idx_;
      }
      if (terms_.find(idx_) == terms_.end()) {
        throw CosaException("Missing term for id " + std::to_string(idx_));
      }

      Term term_ = terms_.at(idx_);
      if (negated_) {
        if (term_->get_sort()->get_sort_kind() == BV) {
          term_ = solver_->make_term(BVNot, term_);
        } else {
          term_ = solver_->make_term(Not, term_);
        }
      }
      termargs_.push_back(term_);
    }

    /******************************** Handle special cases ********************************/
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

      terms_[l_->id] = rts_.make_state(symbol_, linesort_);
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
      terms_[l_->id] = rts_.make_input(symbol_, linesort_);
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
    else if(l_->tag == BTOR2_TAG_sort)
    {
      switch(l_->sort.tag)
      {
      case BTOR2_TAG_SORT_bitvec:
        {
          linesort_=solver_->make_sort(BV, l_->sort.bitvec.width);
          sorts_[l_->id]=linesort_;
          break;
        }
      case BTOR2_TAG_SORT_array:
        {
          linesort_=solver_->make_sort(ARRAY, sorts_.at(l_->sort.array.index), sorts_.at(l_->sort.array.element));
          sorts_[l_->id]=linesort_;
          break;
        }
      default:
        // TODO: maybe only check this in debug? or could always check cause it's really bad
        throw CosaException("Unknown sort tag");
      }
    }
    else if (l_->tag == BTOR2_TAG_constraint)
    {
      rts_.add_constraint(bv_to_bool(termargs_[0]));
    }
    else if (l_->tag == BTOR2_TAG_init)
    {
      rts_.constrain_init(solver_->make_term(Equal, termargs_));
    }
    else if (l_->tag == BTOR2_TAG_next)
    {
      rts_.set_next(termargs_[0], termargs_[1]);
    }
    else if (l_->tag == BTOR2_TAG_bad)
    {
      badvec_.push_back(termargs_[0]);
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
      cval_ = mpz_class(l_->constant, basemap.at(l_->tag));
      terms_[l_->id] = solver_->make_value(cval_.get_str(10), linesort_);
    }
    else if (l_->tag == BTOR2_TAG_one)
    {
      cval_ = mpz_class("1", 10);
      terms_[l_->id] = solver_->make_value(cval_.get_str(10), linesort_);
    }
    else if (l_->tag == BTOR2_TAG_ones)
    {
      cval_ = mpz_class(string(linesort_->get_width(), '1').c_str(), 2);
      terms_[l_->id] = solver_->make_value(cval_.get_str(10), linesort_);
    } else if (l_->tag == BTOR2_TAG_zero) {
      terms_[l_->id] = solver_->make_value(0, linesort_);
    } else if (l_->tag == BTOR2_TAG_slice) {
      terms_[l_->id] = solver_->make_term(Op(Extract, l_->args[1], l_->args[2]), bool_to_bv(termargs_[0]));
    } else if (l_->tag == BTOR2_TAG_sext) {
      terms_[l_->id] = solver_->make_term(Op(Sign_Extend, l_->args[1]), bool_to_bv(termargs_[0]));
    } else if (l_->tag == BTOR2_TAG_uext) {
      terms_[l_->id] = solver_->make_term(Op(Zero_Extend, l_->args[1]), bool_to_bv(termargs_[0]));
    } else if (l_->tag == BTOR2_TAG_rol) {
      terms_[l_->id] = solver_->make_term(Op(Rotate_Left, l_->args[1]), bool_to_bv(termargs_[0]));
    } else if (l_->tag == BTOR2_TAG_ror) {
      terms_[l_->id] = solver_->make_term(Op(Rotate_Right, l_->args[1]), bool_to_bv(termargs_[0]));
    } else if (l_->tag == BTOR2_TAG_inc) {
      Term t = bool_to_bv(termargs_[0]);
      terms_[l_->id] = solver_->make_term(BVAdd, t, solver_->make_value(1, t->get_sort()));
    } else if (l_->tag == BTOR2_TAG_dec) {
      Term t = bool_to_bv(termargs_[0]);
      terms_[l_->id] = solver_->make_term(BVSub, t, solver_->make_value(1, t->get_sort()));
    } else if (l_->tag == BTOR2_TAG_redand) {
      Term t = bool_to_bv(termargs_[0]);
      Term ones = solver_->make_value(pow(2, t->get_sort()->get_width())-1, t->get_sort());
      terms_[l_->id] = solver_->make_term(BVComp, t, ones);
    } else if (l_->tag == BTOR2_TAG_redor) {
      Term t = bool_to_bv(termargs_[0]);
      Term zero = solver_->make_value(0, t->get_sort());
      terms_[l_->id] = solver_->make_term(BVComp, t, zero);
    } else if (l_->tag == BTOR2_TAG_redxor) {
      Term t = bool_to_bv(termargs_[0]);
      unsigned int width = t->get_sort()->get_width();
      Term res = solver_->make_term(Op(Extract, width-1, width-1), t);
      for (int i = width-1; i >= 0; i--)
      {
        res = solver_->make_term(BVXor, res, solver_->make_term(Op(Extract, i, i), t));
      }
      terms_[l_->id] = res;
    } else if (l_->tag == BTOR2_TAG_ite) {
      Term cond = bv_to_bool(termargs_[0]);
      TermVec tv = lazy_convert({termargs_[1], termargs_[2]});
      terms_[l_->id] = solver_->make_term(Ite, cond, tv[0], tv[1]);
    } else if (overflow_ops.find(l_->tag) != overflow_ops.end()) {
      Term t0 = bool_to_bv(termargs_[0]);
      Term t1 = bool_to_bv(termargs_[0]);
      int width = linesort_->get_width();
      t0 = solver_->make_term(Op(Zero_Extend, 1), t0);
      t1 = solver_->make_term(Op(Zero_Extend, 1), t1);
      terms_[l_->id] = solver_->make_term(Op(Extract, width, width),
                                  solver_->make_term(bvopmap.at(l_->tag), t0, t1));
    }
    /******************************** Handle general case
       ********************************/
    else {
      if(boolopmap.find(l_->tag) != boolopmap.end())
      {
        termargs_ = lazy_convert(termargs_);
      }
      else
      {
        for (int i = 0; i < termargs_.size(); i++)
        {
          termargs_[i] = bool_to_bv(termargs_[i]);
        }
      }
      terms_[l_->id] = solver_->make_term(bvopmap.at(l_->tag), termargs_);
    }
  }

  btor2parser_delete(reader_);

}

}
