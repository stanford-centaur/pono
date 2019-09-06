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
  TermVec args;
  unordered_map<int, Sort> sorts;
  unordered_map<int, Term> terms;
  string symbol;

  mpz_class cval;
  const unordered_map<Btor2Tag, int> base(
                                          { {BTOR2_TAG_const, 2},
                                            {BTOR2_TAG_constd, 10},
                                            {BTOR2_TAG_consth, 16} });

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

      fts.name_term(symbol, terms.at(l->id));
    }
    else if (l->constant)
    {
      cval = mpz_class(l->constant, base.at(l->tag));
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


    /******************************** Gather term arguments ********************************/
    // TODO: handle sext etc...
    args.clear();
    args.reserve(l->nargs);
    for(i = 0; i < l->nargs; i++)
    {
      args.push_back(terms.at(l->args[i]));
    }

  }

  btor2parser_delete(reader);

}
