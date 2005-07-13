/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2003 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
 * Tobias Ringström, and Jonas Samuelsson.                                   *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*!
  \file
  \brief Implementation of error handling functions.
  \author Tobias Ringström

  $Revision$

  $Date$
*/

#include "itpp/base/itassert.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdexcept>

using std::cerr;
using std::endl;
using std::flush;
using std::string;

namespace itpp { 

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  static bool warnings_enabled = true;
  static std::ostream *warn = &cerr;
  static bool it_using_exceptions = ITPP_DEFAULT_EXCEPTIONS;

#endif //DOXYGEN_SHOULD_SKIP_THIS

  void it_assert_f(string ass, string msg, string file, int line)
  {
    char line_str[100];

    string error = "*** Assertation failed in ";
    error += file;
    error += " on line ";
    sprintf(line_str, "%d", line);
    error += line_str;
    error += ":\n";
    error += msg;
    error += " (";
    error += ass;
    error += ")";
    cerr << error << endl << flush;
    if (it_using_exceptions)
      throw std::runtime_error(error);
    else
      abort();
  }

  void it_error_f(string msg, string file, int line)
  {
    char line_str[100];

    string error = "*** Error in ";
    error += file;
    error += " on line ";
    sprintf(line_str, "%d", line);
    error += line_str;
    error += ":";
    error += msg;
    cerr << error << endl << flush;
    if (it_using_exceptions)
      throw std::runtime_error(error);
    else
      abort();
  }

  void it_warning_f(string msg, string file, int line)
  {
    if (warnings_enabled)
      (*warn) << "*** Warning in " << file << " on line " << line << ":" << endl
	      << msg << endl << flush;
  }

  void it_enable_exceptions(bool on)
  {
    it_using_exceptions = on;
  }

  void it_enable_warnings()
  {
    warnings_enabled = true;
  }

  void it_disable_warnings()
  {
    warnings_enabled = false;
  }

  void it_redirect_warnings(std::ostream *warn_stream)
  {
    warn = warn_stream;
  }

} //namespace itpp
