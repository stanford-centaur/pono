#pragma once

// use the header only implementation
#define FMT_HEADER_ONLY

#include <fmt/format.h>
#include <iostream>

#include <smt-switch/smt.h>

/****************************** Support for printing smt-switch objects
 * *********************************/

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
    return format_to(ctx.out(), t->to_string());
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
    return format_to(ctx.out(), s->to_string());
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

class Log
{
 public:
  Log() : verbosity(0), verbosity_set(false) {}

  /* Logs to the terminal using Python-style format string
   * @param level the verbosity level to print this log (prints for any
   * verbosity greater than this level)
   * @param format the format string
   * @param args comma separated list of inputs for the format string
   */
  template <typename... Args>
  void log(unsigned int level,
           const std::string & format,
           const Args &... args) const
  {
    if (level <= verbosity)
    {
      std::cout << fmt::format(format, args...) << std::endl;
    }
  }

  /* Logs to the terminal using Python-style format string in a range of
   * verbosities: [lower, upper]
   * @param lower the lower cutoff for verbosity
   * @param upper the upper cutoff for verbosity
   * @param format the format string
   * @param args comma separated list of inputs for the format string
   */
  template <typename... Args>
  void log(unsigned int lower,
           unsigned int upper,
           const std::string & format,
           const Args &... args) const
  {
    if ((lower <= verbosity) && (verbosity <= upper))
    {
      std::cout << fmt::format(format, args...) << std::endl;
    }
  }

  /* set verbosity -- can only be set once
   * @param v the verbosity to set
   */
  void set_verbosity(unsigned int v)
  {
    if (!verbosity_set)
    {
      verbosity = v;
    }
    else
    {
      throw CosaException("Can only set logger verbosity once.");
    }
  }

 protected:
  int verbosity;
  bool verbosity_set;
};

Log logger;
