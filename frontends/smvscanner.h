#pragma once

#include <stdio.h>
#include <iostream>

#ifndef yyFlexLexerOnce
#include <FlexLexer.h>
#endif

#include "smvparser.h"

namespace pono {
class SMVEncoder;
class SMVscanner : public yyFlexLexer
{
 public:
  SMVscanner(SMVEncoder & encoder) : _encoder(encoder) {}
  virtual ~SMVscanner() {}
  virtual pono::smvparser::symbol_type yylex(SMVEncoder & encoder);

 private:
  SMVEncoder & _encoder;
};
}  // namespace pono
