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

#pragma once

#include <iostream>
#include <utility>

// pono's CMake target already requests the header-only build of fmt, but
// define it here too so that including the installed copy of this header needs
// nothing but fmt's headers. Guarded because the target defines it with a
// value.
#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY 1
#endif

#include <fmt/format.h>

#include "smt-switch/smt.h"
#include "utils/exceptions.h"

/****************************** Support for printing smt-switch objects
 * *********************************/

/** Formats an smt-switch object by delegating to the std::string formatter.
 *  Delegating means the printed text is passed as an argument rather than as a
 *  format string, so terms containing curly brackets print verbatim, and the
 *  inherited parse() accepts the usual string format specs.
 *  @param type the smt-switch type to format
 *  @param expr an expression producing the text for a value named v
 */
#define PONO_STRING_FORMATTER(type, expr)                   \
  template <>                                               \
  struct fmt::formatter<type> : formatter<std::string>      \
  {                                                         \
    auto format(const type & v, format_context & ctx) const \
    {                                                       \
      return formatter<std::string>::format(expr, ctx);     \
    }                                                       \
  };

PONO_STRING_FORMATTER(smt::Term, v->to_string())
PONO_STRING_FORMATTER(smt::Sort, v->to_string())
PONO_STRING_FORMATTER(smt::PrimOp, smt::to_string(v))
PONO_STRING_FORMATTER(smt::Op, v.to_string())
PONO_STRING_FORMATTER(smt::Result, v.to_string())

/*********************** End overloaded methods for printing smt-switch objects
 * **********************/

/*************************************** Logger class
 * ************************************************/
// Meant to be used as a singleton class -- instantiated as logger below

namespace pono {

class Log
{
 public:
  Log() : verbosity(0), verbosity_set(false) {}

  Log(size_t v) : verbosity(v), verbosity_set(true) {}

  /* Logs to the output stream using Python-style format string
   * @param level the verbosity level to print this log (prints for any
   * verbosity greater than this level)
   * @param format the format string
   * @param args comma separated list of inputs for the format string
   */
  template <typename... Args>
  void log_to_stream(size_t level,
                     std::ostream & output_stream,
                     fmt::format_string<Args...> format,
                     Args &&... args) const
  {
    if (level <= verbosity) {
      output_stream << fmt::format(format, std::forward<Args>(args)...)
                    << std::endl;
    }
  }

  /* Logs to the terminal using Python-style format string
   * @param level the verbosity level to print this log (prints for any
   * verbosity greater than this level)
   * @param format the format string
   * @param args comma separated list of inputs for the format string
   */
  template <typename... Args>
  void log(size_t level,
           fmt::format_string<Args...> format,
           Args &&... args) const
  {
    // Formatted here rather than by delegating to log_to_stream, because
    // passing a format_string on to a second checked overload would deduce a
    // different Args pack and reject every call site.
    if (level <= verbosity) {
      std::cerr << fmt::format(format, std::forward<Args>(args)...)
                << std::endl;
    }
  }

  /* set verbosity -- can only be set once
   * @param v the verbosity to set
   */
  void set_verbosity(size_t v)
  {
    if (!verbosity_set) {
      verbosity = v;
    } else {
      throw PonoException("Can only set logger verbosity once.");
    }
  }

 protected:
  size_t verbosity;
  bool verbosity_set;
};

// globally available logger instance
extern Log logger;

void set_global_logger_verbosity(size_t v);

}  // namespace pono
