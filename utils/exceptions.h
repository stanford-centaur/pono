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


#include <exception>
#include <string>

/**
   Base exception for this codebase.

   All other exceptions used in the code should be derived
     classes of this.
 */

class PonoException : public std::exception
{
 public:
  /** Constructor (C strings).
   *  @param message C-style string error message.
   *                 The string contents are copied upon construction.
   *                 Hence, responsibility for deleting the char* lies
   *                 with the caller.
   */
  explicit PonoException(const char * message) : msg(message) {}

  /** Constructor (C++ STL strings).
   *  @param message The error message.
   */
  explicit PonoException(const std::string & message) : msg(message) {}

  /** Destructor.
   * Virtual to allow for subclassing.
   */
  virtual ~PonoException() throw() {}

  /** Returns a pointer to the (constant) error description.
   *  @return A pointer to a const char*. The underlying memory
   *          is in posession of the Exception object. Callers must
   *          not attempt to free the memory.
   */
  virtual const char * what() const throw() { return msg.c_str(); }

 protected:
  /** Error message.
   */
  std::string msg;
};

