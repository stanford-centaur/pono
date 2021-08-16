/*********************                                                        */
/*! \file vmt_encoder.h
** \verbatim
** Top contributors (to current version):
**   Makai Mann
** This file is part of the pono project.
** Copyright (c) 2019 by the authors listed in the file AUTHORS
** in the top-level source directory) and their institutional affiliations.
** All rights reserved.  See the file LICENSE in the top-level source
** directory for licensing information.\endverbatim
**
** \brief Frontend for the Verification Modulo Theories (VMT) format.
**        See https://vmt-lib.fbk.eu/ for more information.
**
**
**/

#pragma once

#include <iostream>

#include "assert.h"
#include "core/rts.h"
#include "smt-switch/smt.h"
#include "smt-switch/smtlib_reader.h"
#include "utils/exceptions.h"

namespace pono {
class VMTEncoder : public smt::SmtLibReader
{
 public:
  VMTEncoder(std::string filename, RelationalTransitionSystem & rts);

  typedef SmtLibReader super;

  void new_symbol(const std::string & name, const smt::Sort & sort) override;

  void term_attribute(const smt::Term & term,
                      const std::string & keyword,
                      const std::string & value) override;

  const smt::TermVec & propvec() const { return propvec_; }

 protected:
  std::string filename_;

  RelationalTransitionSystem & rts_;

  smt::TermVec propvec_;
};

}  // namespace pono
