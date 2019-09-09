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

#include "exceptions.h"
#include "fts.h"
#include "smt-switch/smt.h"

namespace cosa
{
  class BTOR2Encoder
  {
  public:
    BTOR2Encoder(std::string filename, FunctionalTransitionSystem & fts)
      : fts_(fts),
        solver_(fts.solver())
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
    smt::SmtSolver & solver_;
    cosa::FunctionalTransitionSystem & fts_;

    // Useful variables
    smt::Sort linesort_;
    smt::TermVec termargs_;
    std::unordered_map<int, smt::Sort> sorts_;
    std::unordered_map<int, smt::Term> terms_;
    std::string symbol_;
    mpz_class cval_;

    smt::TermVec badvec_;
    smt::TermVec justicevec_;
    smt::TermVec fairvec_;

    Btor2Parser* reader_;
    Btor2LineIterator it_;
    Btor2Line* l_;
    size_t i_;

  };
}

#endif
