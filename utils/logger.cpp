/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Makai Mann
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

#include "utils/logger.h"

std::string remove_curly_brackets(std::string s)
{
  std::size_t pos;
  while ((pos = s.find("{")) != std::string::npos) {
    s.replace(pos, 1, "");
  }
  while ((pos = s.find("}")) != std::string::npos) {
    s.replace(pos, 1, "");
  }
  return s;
}

// declare a global logger
namespace pono {
Log logger;

void set_global_logger_verbosity(size_t v) { logger.set_verbosity(v); }
}  // namespace pono
