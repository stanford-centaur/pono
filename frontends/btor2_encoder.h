#ifndef BTOR2_ENCODER_H
#define BTOR2_ENCODER_H

extern "C"
{
#include "btor2parser/btor2parser.h"
}
#include "gmpxx.h"

#include "assert.h"
#include <iostream>
#include "math.h"
#include <stdio.h>
#include <string>
#include <unordered_map>

#include "fts.h"
#include "smt-switch/smt.h"

namespace cosa
{
  class BTOR2Encoder
  {
  public:
    BTOR2Encoder(std::string filename, FunctionalTransitionSystem & fts)
      : fts(fts),
        s(fts.solver())
     {
       parse(filename);
     };

  protected:

    // converts booleans to bitvector of size one
    smt::Term bool_to_bv(smt::Term t);
    // converts bitvector of size one to boolean
    smt::Term bv_to_bool(smt::Term t);
    // lazy conversion
    // takes a list of booleans / bitvectors of size one
    // and lazily converts them to the majority
    smt::TermVec lazy_convert(smt::TermVec);

    // parse a btor2 file
    void parse(std::string filename);

    // Important members
    smt::SmtSolver & s;
    cosa::FunctionalTransitionSystem & fts;

    // Useful variables
    smt::Sort linesort;
    smt::TermVec termargs;
    std::unordered_map<int, smt::Sort> sorts;
    std::unordered_map<int, smt::Term> terms;
    std::string symbol;
    mpz_class cval;

    smt::TermVec badvec;
    smt::TermVec justicevec;
    smt::TermVec fairvec;

    Btor2Parser* reader;
    Btor2LineIterator it;
    Btor2Line* l;
    size_t i;

  };
}

#endif
