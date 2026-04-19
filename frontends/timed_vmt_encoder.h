#pragma once

#include <iostream>

#include "assert.h"
#include "core/rts.h"
#include "smt-switch/smt.h"
#include "smt-switch/smtlib_reader.h"
#include "utils/exceptions.h"
#include "vmt_encoder.h"
#include "core/tts.h"

namespace pono {
class TimedVMTEncoder : public smt::SmtLibReader
{
 public:
  TimedVMTEncoder(std::string filename, TimedTransitionSystem & tts);

  void new_symbol(const std::string & name, const smt::Sort & sort) override;

  void term_attribute(const smt::Term & term,
                      const std::string & keyword,
                      const std::string & value) override;
  const smt::TermVec & propvec() const { return propvec_; }

  protected:
  TimedTransitionSystem & tts_;
  std::string filename_;
  smt::TermVec propvec_;  
};

}  // namespace pono
