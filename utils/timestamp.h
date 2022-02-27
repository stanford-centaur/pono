/*********************                                                        */
/*! \file
 ** \verbatim
 ** Top contributors (to current version):
 **   Florian Lonsing
 ** This file is part of the pono project.
 ** Copyright (c) 2021 by the authors listed in the file AUTHORS
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
#include <chrono>

/*************************************** time stamp functions
 * ************************************************/

namespace pono {

// Inspiration:
// https://levelup.gitconnected.com/8-ways-to-measure-execution-time-in-c-c-48634458d0f9
//
// Note: frequently taking time stamps, e.g. in a loop, may be costly.
// We measure timestamp in nanoseconds and then convert, e.g., to seconds.
  
typedef std::chrono::system_clock pono_clock;
typedef std::chrono::time_point<pono_clock> pono_time_stamp;
typedef std::chrono::duration<long long int, std::nano> pono_time_duration;

// take current time stamp
static pono_time_stamp timestamp () {
  return pono_clock::now();
}

// compute duration in nanoseconds between two given time stamps
static pono_time_duration timestamp_diff (
  pono_time_stamp begin, pono_time_stamp end) {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
}

// convert duration in nanoseconds computed by 'timestamp_diff' to a string
static string time_duration_to_sec_string (pono_time_duration d) {
  ostringstream out;
  out << d.count() * 1e-9;
  return out.str();
}

}  // namespace pono
