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
 ** \brief some shortcut for STL containers
 **
 ** 
 **/

#pragma once

#include <algorithm>
#include <set>
#include <string>
#include <unordered_set>

// actually no use, because micros are not bound by namespace
// should only be included in .cpp
namespace cosa {

#define UNION(a, b, r)                                                         \
  (std::set_union((a).begin(), (a).end(), (b).begin(), (b).end(),              \
                  std::inserter((r), (r).end())))
#define INTERSECT(a, b, r)                                                     \
  (std::set_intersection((a).begin(), (a).end(), (b).begin(), (b).end(),       \
                         std::inserter((r), (r).end())))
#define DIFFERENCE(a, b, r)                                                    \
  (std::set_difference((a).begin(), (a).end(), (b).begin(), (b).end(),         \
                       std::inserter((r), (r).end())))
#define SYMDIFF(a, b, r)                                                       \
  (std::set_symmetric_difference((a).begin(), (a).end(), (b).begin(),          \
                                 (b).end(), std::inserter((r), (r).end())))

#define IN(e, s) ((s).find(e) != (s).end())
#define IN_p(e, s) ((s)->find(e) != (s)->end())

#define S_IN(sub, s) ((s).find(sub) != (s).npos)

#define FIND_IN(e,s) ((std::find((s).begin(), (s).end(), (e))) != (s).end())

template<typename MAP>
const typename MAP::mapped_type& get_with_default(const MAP& m, 
                                             const typename MAP::key_type& key, 
                                             const typename MAP::mapped_type& defval)
{
    typename MAP::const_iterator it = m.find(key);
    if (it == m.end())
        return defval;

    return it->second;
}

template<class T> bool is_union_empty(
    const std::unordered_set<T> & a, 
    const std::unordered_set<T> & b) {

  const std::unordered_set<T> & small = a.size() < b.size() ? a : b;
  const std::unordered_set<T> & big = a.size() < b.size() ? b : a;
  for(auto && p : small) {
    if (big.find(p) != big.end())
      return false;
  }
  return true;
}

}  // namespace cosa