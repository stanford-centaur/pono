/*********************                                                        */
/*! \file 
 ** \verbatim
 ** Top contributors (to current version):
 **   Hongce Zhang
 ** This file is part of the cosa2 project.
 ** Copyright (c) 2019 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief 
 **
 ** 
 **/

#include "utils/logger.h"
#include "frontends/btor2_encoder.h"
#include "smt-switch/boolector_factory.h"

#include "vcd_witness_printer.h"

#include <iostream>

using namespace smt;
using namespace std;

namespace cosa {

const std::string_view vcd_header (R"**(
$date %date% $end
$version COSA2 $end
$timescale 1 ns $end
)**");

VCDWitnessPrinter::VCDWitnessPrinter(const BTOR2Encoder & btor_enc,
                    const std::vector<smt::UnorderedTermMap> & cex) :
  inputs(btor_enc.inputsvec()),
  states(btor_enc.statesvec()),
  no_next_states(btor_enc.no_next_states()),
  has_states_without_next(!no_next_states.empty())
{
  // figure out the variables and their scopes

  // you also need to figure out the indices used for each array
  // and dump these values when needed and dump the default value
  logger.log(0, "-------------Input Dump-----------------");
  for (auto && i: inputs) {
    logger.log(0, "{}", i->to_string());
  }
  logger.log(0, "-------------State Dump-----------------");
  for (auto && s: states) {
    logger.log(0, "{}", s->to_string());
  }
  logger.log(0, "-------------No Next States-----------------");
  for (auto && s: no_next_states) {
    logger.log(0, "{}", s.second->to_string());
  }
  for (uint64_t fidx = 0; fidx < cex.size(); ++ fidx) {
    logger.log(0, "------------- CEX : F{} -----------------", fidx);
    for (auto && t : cex.at(fidx)) {
      logger.log(0, "{} -> {}", t.first->to_string(), t.second->to_string() );
    }
  }
}



} // namespace cosa
