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
    Sort bv1sort = s->make_sort(BV, 1);
    return s->make_term(Ite, s->make_value(1, bv1sort), s->make_value(0, bv1sort));
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
      throw "Can't convert non-width 1 bitvector to bool.";
    }
    return s->make_term(Equal, t, s->make_value(1, s->make_sort(BV, 1)));
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
    throw "Could not open " + filename;
  }

  reader = btor2parser_new();

  if(!btor2parser_read_lines(reader, input_file))
  {
    throw "Error parsing btor file.";
  }

  it = btor2parser_iter_init(reader);
  while((l = btor2parser_iter_next(&it)))
  {
    /******************************** Identify sort ********************************/
    if (l->tag != BTOR2_TAG_sort && l->sort.id)
    {
      linesort=sorts.at(l->sort.id);
    }

    /******************************** Gather term arguments ********************************/
    termargs.clear();
    termargs.reserve(l->nargs);
    for(i = 0; i < l->nargs; i++)
    {
      if (terms.find(l->args[i]) == terms.end())
      {
        break;
      }
      termargs.push_back(terms.at(l->args[i]));
    }

    /******************************** Handle special cases ********************************/
    if (l->tag == BTOR2_TAG_state)
    {
      if (l->symbol)
      {
        symbol = l->symbol;
      }
      else
      {
        symbol = "state" + to_string(l->id);
      }

      terms[l->id] = fts.make_state(symbol, linesort);
    }
    else if (l->tag == BTOR2_TAG_input)
    {
      if (l->symbol)
      {
        symbol = l->symbol;
      }
      else
      {
        symbol = "input" + to_string(l->id);
      }

      terms[l->id] = fts.make_input(symbol, linesort);
    }
    else if (l->tag == BTOR2_TAG_output)
    {
      if (l->symbol)
      {
        symbol = l->symbol;
      }
      else
      {
        symbol = "output" + to_string(l->id);
      }

      fts.name_term(symbol, termargs[0]);
    }
    else if(l->tag == BTOR2_TAG_sort)
    {
      switch(l->sort.tag)
      {
      case BTOR2_TAG_SORT_bitvec:
        {
          linesort=s->make_sort(BV, l->sort.bitvec.width);
          sorts[l->id]=linesort;
          break;
        }
      case BTOR2_TAG_SORT_array:
        {
          linesort=s->make_sort(ARRAY, sorts.at(l->sort.array.index), sorts.at(l->sort.array.element));
          sorts[l->id]=linesort;
          break;
        }
      default:
        // TODO: maybe only check this in debug? or could always check cause it's really bad
        throw "Unknown sort tag";
      }
    }
    else if (l->tag == BTOR2_TAG_constraint)
    {
      fts.add_constraint(bv_to_bool(termargs[0]));
    }
    else if (l->tag == BTOR2_TAG_init)
    {
      fts.constrain_init(s->make_term(Equal, termargs));
    }
    else if (l->tag == BTOR2_TAG_next)
    {
      fts.set_next(termargs[0], termargs[1]);
    }
    else if (l->tag == BTOR2_TAG_bad)
    {
      badvec.push_back(termargs[0]);
    }
    else if (l->tag == BTOR2_TAG_justice)
    {
      std::cout << "Warning: ignoring justice term" << std::endl;
      justicevec.push_back(termargs[0]);
    }
    else if (l->tag == BTOR2_TAG_fair)
    {
      std::cout << "Warning: ignoring fair term" << std::endl;
      fairvec.push_back(termargs[0]);
    }
    else if (l->constant)
    {
      cval = mpz_class(l->constant, basemap.at(l->tag));
      terms[l->id] = s->make_value(cval.get_str(10), linesort);
    }
    else if (l->tag == BTOR2_TAG_one)
    {
      cval = mpz_class("1", 10);
      terms[l->id] = s->make_value(cval.get_str(10), linesort);
    }
    else if (l->tag == BTOR2_TAG_ones)
    {
      cval = mpz_class(string(linesort->get_width(), '1').c_str(), 2);
      terms[l->id] = s->make_value(cval.get_str(10), linesort);
    }
    else if (l->tag == BTOR2_TAG_slice)
    {
      terms[l->id] = s->make_term(Op(Extract, l->args[1], l->args[2]), bool_to_bv(termargs[0]));
    }
    else if (l->tag == BTOR2_TAG_sext)
    {
      terms[l->id] = s->make_term(Op(Sign_Extend, l->args[1]), bool_to_bv(termargs[0]));
    }
    else if (l->tag == BTOR2_TAG_uext)
    {
      terms[l->id] = s->make_term(Op(Zero_Extend, l->args[1]), bool_to_bv(termargs[0]));
    }
    else if (l->tag == BTOR2_TAG_rol)
    {
      terms[l->id] = s->make_term(Op(Rotate_Left, l->args[1]), bool_to_bv(termargs[0]));
    }
    else if (l->tag == BTOR2_TAG_ror)
    {
      terms[l->id] = s->make_term(Op(Rotate_Right, l->args[1]), bool_to_bv(termargs[0]));
    }
    else if (l->tag == BTOR2_TAG_inc)
    {
      Term t = bool_to_bv(termargs[0]);
      terms[l->id] = s->make_term(BVAdd, t, s->make_value(1, t->get_sort()));
    }
    else if (l->tag == BTOR2_TAG_dec)
    {
      Term t = bool_to_bv(termargs[0]);
      terms[l->id] = s->make_term(BVSub, t, s->make_value(1, t->get_sort()));
    }
    else if (l->tag == BTOR2_TAG_redand)
    {
      Term t = bool_to_bv(termargs[0]);
      Term ones = s->make_value(pow(2, t->get_sort()->get_width())-1, t->get_sort());
      terms[l->id] = s->make_term(BVComp, t, ones);
    }
    else if (l->tag == BTOR2_TAG_redor)
    {
      Term t = bool_to_bv(termargs[0]);
      Term zero = s->make_value(0, t->get_sort());
      terms[l->id] = s->make_term(BVComp, t, zero);
    }
    else if (l->tag == BTOR2_TAG_redxor)
    {
      Term t = bool_to_bv(termargs[0]);
      unsigned int width = t->get_sort()->get_width();
      Term res = s->make_term(Op(Extract, width-1, width-1), t);
      for (int i = width-1; i >= 0; i--)
      {
        res = s->make_term(BVXor, res, s->make_term(Op(Extract, i, i), t));
      }
      terms[l->id] = res;
    }
    else if (l->tag == BTOR2_TAG_ite)
    {
      Term cond = bv_to_bool(termargs[0]);
      TermVec tv = lazy_convert({termargs[1], termargs[2]});
      terms[l->id] = s->make_term(Ite, cond, tv[0], tv[1]);
    }
    else if (overflow_ops.find(l->tag) != overflow_ops.end())
    {
      Term t0 = bool_to_bv(termargs[0]);
      Term t1 = bool_to_bv(termargs[0]);
      int width = linesort->get_width();
      t0 = s->make_term(Op(Zero_Extend, 1), t0);
      t1 = s->make_term(Op(Zero_Extend, 1), t1);
      terms[l->id] = s->make_term(Op(Extract, width, width),
                                  s->make_term(bvopmap.at(l->tag), t0, t1));
    }
    /******************************** Handle general case ********************************/
    else
    {
      if(boolopmap.find(l->tag) != boolopmap.end())
      {
        termargs = lazy_convert(termargs);
      }
      else
      {
        for (int i = 0; i < termargs.size(); i++)
        {
          termargs[i] = bool_to_bv(termargs[i]);
        }
      }
      terms[l->id] = s->make_term(bvopmap.at(l->tag), termargs);
    }

  }

  btor2parser_delete(reader);

}

}
