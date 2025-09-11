/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann, Ahmed Irfan, Po-Chun Chien, Áron Ricardo Perez-Lopez
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

#include <iostream>
#include <string>
#include <vector>

#include "core/ts.h"
#include "frontends/btor2_encoder.h"
#include "smt-switch/smt.h"

namespace pono {
void print_witness_btor(const BTOR2Encoder & btor_enc,
                        const std::vector<smt::UnorderedTermMap> & cex,
                        const TransitionSystem & ts,
                        std::ostream & output_stream = std::cout);

void dump_witness_btor(const BTOR2Encoder & btor_enc,
                       std::vector<smt::UnorderedTermMap> & cex,
                       const TransitionSystem & ts,
                       const unsigned int prop_idx,
                       const std::string & witness_filename);
}  // namespace pono
