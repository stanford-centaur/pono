/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Hongce Zhang
 ** This file is part of the pono project.
 ** Copyright (c) 2020 by the authors listed in the file AUTHORS
 ** in the top-level source directory) and their institutional affiliations.
 ** All rights reserved.  See the file LICENSE in the top-level source
 ** directory for licensing information.\endverbatim
 **
 ** \brief SyGuS string helper
 **
 **
 **/

#include "str_util.h"

#include <cassert>
#include <iostream>
#include <sstream>

namespace pono {

namespace syntax_analysis {

// it is of course possible to update it with arbitrary base
std::string IntToStrCustomBase(uint64_t value, unsigned base, bool uppercase)
{
  assert(base > 1 && base <= 36);
  if (value == 0) return "0";
  std::string ret;
  while (value != 0) {
    unsigned digit_val = value % base;
    char digit = (digit_val < 10) ? ('0' + digit_val)
                                  : ((uppercase ? 'A' : 'a') + digit_val - 10);
    ret = digit + ret;
    value /= base;
  }
  return ret;
}

unsigned long long StrToULongLong(const std::string & str, int base)
{
  return std::stoull(str, NULL, base);
}

/// Trim a string from start (in place)
void StrLeftTrim(std::string & s)
{
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
#if defined(_WIN32) || defined(_WIN64)
            return !std::isspace(ch, std::locale("en_US.UTF8"));
#else
        return !std::isspace(static_cast<unsigned char>(ch));
#endif
          }));
}

/// Trim a string from end (in place)
void StrRightTrim(std::string & s)
{
  s.erase(std::find_if(s.rbegin(),
                       s.rend(),
                       [](int ch) {
#if defined(_WIN32) || defined(_WIN64)
                         return !std::isspace(ch, std::locale("en_US.UTF8"));
#else
        return !std::isspace(static_cast<unsigned char>(ch));
#endif
                       })
              .base(),
          s.end());
}

/// Trim a string from both ends (in place)
void StrTrim(std::string & s)
{
  StrLeftTrim(s);
  StrRightTrim(s);
}

std::vector<std::string> Split(const std::string & str,
                               const std::string & delim)
{
  std::vector<std::string> tokens;
  size_t prev = 0, pos = 0;
  do {
    pos = str.find(delim, prev);
    if (pos == std::string::npos) pos = str.length();
    std::string token = str.substr(prev, pos - prev);
    if (!token.empty()) tokens.push_back(token);
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());
  return tokens;
}

std::vector<std::string> SplitSpaceTabEnter(const std::string & str)
{
  std::vector<std::string> result;
  std::istringstream iss(str);
  for (std::string s; iss >> s;) result.push_back(s);
  return result;
}

std::string Join(const std::vector<std::string> & in, const std::string & delim)
{
  std::string ret;
  std::string d = "";
  for (auto && s : in) {
    ret += (d + s);
    d = delim;
  }
  return ret;
}

/// Remove whitespace ' \n\t\r\f\v'
std::string RemoveWhiteSpace(const std::string & in)
{
  auto s = in;
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  return s;
}

/// Replace all occurrences of substring a by substring b
std::string ReplaceAll(const std::string & str,
                       const std::string & a,
                       const std::string & b)
{
  std::string result;
  size_t find_len = a.size();
  size_t pos, from = 0;
  while (std::string::npos != (pos = str.find(a, from))) {
    result.append(str, from, pos - from);
    result.append(b);
    from = pos + find_len;
  }
  result.append(str, from, std::string::npos);
  return result;
}

bool StrEndsWith(const std::string & str, const std::string & suffix)
{
  return str.size() >= suffix.size()
         && 0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

bool StrStartsWith(const std::string & str, const std::string & prefix)
{
  return str.size() >= prefix.size()
         && 0 == str.compare(0, prefix.size(), prefix);
}

void mul2(std::vector<char> & v)
{
  char carry = 0;
  for (auto pos = v.begin(); pos != v.end(); ++pos) {
    *pos = (*pos) * 2 + carry;
    if (*pos >= 10) {
      carry = *pos / 10;
      *pos = *pos % 10;
    } else
      carry = 0;
  }
  if (carry) v.push_back(carry);
}

void add1(std::vector<char> & v)
{
  char carry = 1;
  for (auto pos = v.begin(); pos != v.end(); ++pos) {
    *pos = *pos + carry;
    if (*pos >= 10) {
      carry = *pos / 10;
      *pos = *pos % 10;
    } else
      carry = 0;
  }
  if (carry) v.push_back(carry);
}

}  // namespace syntax_analysis

const std::string & lookup_or_key(
    const std::unordered_map<std::string, std::string> & map,
    const std::string & key)
{
  auto it = map.find(key);
  return (it != map.end() && !it->second.empty()) ? it->second : key;
}

}  // namespace pono
