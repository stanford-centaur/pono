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


#include <iterator>
#include <map>
#include <sstream>
#include <vector>

#include "gmpxx.h"

#include "smt-switch/smt.h"

#include "utils/logger.h"

namespace cosa {

struct VCDScope {
  std::vector<VCDScope> subscopes;
  smt::TermVec variables;
}; // struct VCDScope

class VCDWitnessPrinter {
public:
  // types
  typedef std::map<std::string, std::vector<std::string>> per_mem_indices;
protected:
  const smt::TermVec inputs;
  const smt::TermVec states;
  const std::map<uint64_t, smt::Term> no_next_states;
  const bool has_states_without_next;
  
  VCDScope root_scope;

  // given a name like a.b.c, find the right scope and
  // create if it does not exists
  void check_insert_scope(const std::string& full_name);

public:
  VCDWitnessPrinter(const BTOR2Encoder & btor_enc,
                    const std::vector<smt::UnorderedTermMap> & cex);

}; // class VCDWitnessPrinter

}  // namespace cosa