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

#include <algorithm>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace pono {

namespace syntax_analysis {

/// Trim a string from start (in place)
void StrLeftTrim(std::string & s);

/// Trim a string from end (in place)
void StrRightTrim(std::string & s);

/// Trim a string from both ends (in place)
void StrTrim(std::string & s);

/// Python-style split, return a vector of split strings
std::vector<std::string> Split(const std::string & str,
                               const std::string & delim);

/// Python-style split behavior, delim: space tab enter and their combinations
std::vector<std::string> SplitSpaceTabEnter(const std::string & str);

// (itoa is not part of the standard actually,
// so I hesitated whether to use it actually)
// on the other hand, snprintf only supports 8/10/16
/// Transform int to string with different bases
std::string IntToStrCustomBase(uint64_t value, unsigned base, bool uppercase);

/// Return the value represented in the string in unsigned long long, e.g. "10".
unsigned long long StrToULongLong(const std::string & str, int base);

/// Python-style join, return a string that joins the list by the delim
std::string Join(const std::vector<std::string> & in,
                 const std::string & delim);

/// Remove whitespace " \n\t\r\f\v" from the input string
std::string RemoveWhiteSpace(const std::string & in);

/// Replace all occurrences of substring a by substring b
std::string ReplaceAll(const std::string & str,
                       const std::string & a,
                       const std::string & b);

/// Finds out if str ends with suffix
bool StrEndsWith(const std::string & str, const std::string & suffix);

/// Finds out if str starts with prefix
bool StrStartsWith(const std::string & str, const std::string & prefix);

/// multi-precision number : multiply by 2
void mul2(std::vector<char> & v);

/// multi-precision number : add by 1
void add1(std::vector<char> & v);

}  // namespace syntax_analysis

// loop up the key in the map
// if found and value is non-empty, return the value
// otherwise, return the key
const std::string & lookup_or_key(
    const std::unordered_map<std::string, std::string> & map,
    const std::string & key);

}  // namespace pono
