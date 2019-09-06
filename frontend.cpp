extern "C"
{
  #include "btor2parser/btor2parser.h"
}
#include "gmpxx.h"

#include "assert.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include <unordered_map>

#include "fts.h"
#include "smt-switch/smt.h"
#include "smt-switch/boolector_factory.h"

using namespace std;
using namespace smt;
using namespace cosa;

int main(int argc, char** argv)
{
  SmtSolver s = BoolectorSolverFactory::create();
  s->set_opt("produce-models", true);
  FunctionalTransitionSystem fts(s);

  Sort linesort;
  TermVec termargs;
  unordered_map<int, Sort> sorts;
  unordered_map<int, Term> terms;
  string symbol;

  mpz_class cval;
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
  const unordered_map<Btor2Tag, PrimOp> opmap(
                                              {{},
                                              }
    );

  TermVec badvec;
  TermVec justicevec;

  Btor2Parser* reader;
  Btor2LineIterator it;
  Btor2Line* l;
  size_t i;

  if (argc != 2)
  {
    throw "Incorrect number of arguments, takes one btor file.";
  }

  FILE* input_file = fopen(argv[1], "r");

  reader = btor2parser_new();

  if(!btor2parser_read_lines(reader, input_file))
  {
    throw "Error parsing btor file.";
  }

  it = btor2parser_iter_init(reader);
  while((l = btor2parser_iter_next(&it)))
  {
    /******************************** Identify sort ********************************/
    if(l->tag == BTOR2_TAG_sort)
    {
      // TODO: handle the sort tags
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
    else if (l->sort.id)
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
      justicevec.push_back(termargs[0]);
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
      terms[l->id] = s->make_term(Op(Extract, l->args[1], l->args[2]), termargs[0]);
    }
    else if (l->tag == BTOR2_TAG_sext)
    {
      terms[l->id] = s->make_term(Op(Sign_Extend, l->args[1]), termargs[0]);
    }
    else if (overflow_ops.find(l->tag) != overflow_ops.end())
    {
      // TODO: finish this
      int width = linesort->get_width();
      Term t0 = s->make_term(Op(Zero_Extend, 1), termargs[0]);
      Term t1 = s->make_term(Op(Zero_Extend, 1), termargs[1]);
      terms[l->id] = s->make_term(Op(Extract, width, width),
                                  s->make_term(opmap.at(l->tag), t0, t1));
    }
    /******************************** Handle general case ********************************/
    else
    {
      terms[l->id] = s->make_term(opmap.at(l->tag), termargs);
    }

  }

  btor2parser_delete(reader);

}
