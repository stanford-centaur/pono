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

// use the header only implementation
#define FMT_HEADER_ONLY

#include <iostream>

#include "fmt/format.h"
#include "smt-switch/smt.h"
#include "utils/exceptions.h"

/****************************** Support for printing smt-switch objects
 * *********************************/

/** Takes a string and removes the curly brackets
 *  Since fmt/format.h uses {} to denote an argument
 *  it breaks if there are curly brackets in the string
 *  @param s the string to start with
 *  @return the same string but without curly brackets
 */
std::string remove_curly_brackets(std::string s);

// Term
template <>
struct fmt::formatter<smt::Term>
{
  template <typename ParseContext>
  constexpr auto parse(ParseContext & ctx)
  {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(const smt::Term & t, FormatContext & ctx)
  {
    return format_to(ctx.out(), remove_curly_brackets(t->to_string()));
  }
};

// Sort
template <>
struct fmt::formatter<smt::Sort>
{
  template <typename ParseContext>
  constexpr auto parse(ParseContext & ctx)
  {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(const smt::Sort & s, FormatContext & ctx)
  {
    return format_to(ctx.out(), remove_curly_brackets(s->to_string()));
  }
};

// PrimOp
template <>
struct fmt::formatter<smt::PrimOp>
{
  template <typename ParseContext>
  constexpr auto parse(ParseContext & ctx)
  {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(const smt::PrimOp & po, FormatContext & ctx)
  {
    return format_to(ctx.out(), smt::to_string(po));
  }
};

// Op
template <>
struct fmt::formatter<smt::Op>
{
  template <typename ParseContext>
  constexpr auto parse(ParseContext & ctx)
  {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(const smt::Op & o, FormatContext & ctx)
  {
    return format_to(ctx.out(), o.to_string());
  }
};

// Result
template <>
struct fmt::formatter<smt::Result>
{
  template <typename ParseContext>
  constexpr auto parse(ParseContext & ctx)
  {
    return ctx.begin();
  }

  template <typename FormatContext>
  auto format(const smt::Result & r, FormatContext & ctx)
  {
    return format_to(ctx.out(), r.to_string());
  }
};

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
                     const std::string & format,
                     const Args &... args) const
  {
    if (level <= verbosity) {
      output_stream << fmt::format(format, args...) << std::endl;
    }
  }

  /* Logs to the terminal using Python-style format string
   * @param level the verbosity level to print this log (prints for any
   * verbosity greater than this level)
   * @param format the format string
   * @param args comma separated list of inputs for the format string
   */
  template <typename... Args>
  void log(size_t level, const std::string & format, const Args &... args) const
  {
    log_to_stream(level, std::cerr, format, args...);
  }

  /* Logs to the terminal using Python-style format string in a range of
   * verbosities: [lower, upper]
   * @param lower the lower cutoff for verbosity
   * @param upper the upper cutoff for verbosity
   * @param format the format string
   * @param args comma separated list of inputs for the format string
   */
  template <typename... Args>
  void log(size_t lower,
           size_t upper,
           const std::string & format,
           const Args &... args) const
  {
    if ((lower <= verbosity) && (verbosity <= upper)) {
      std::cerr << fmt::format(format, args...) << std::endl;
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
