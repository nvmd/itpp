/*!
 * \file 
 * \brief Implementation of error handling functions
 * \author Tobias Ringstrom
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * -------------------------------------------------------------------------
 */

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <itpp/base/itassert.h>
#include <iostream>
#include <string>
#include <stdexcept>


namespace itpp {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  static bool warnings_enabled = true;
  static std::ostream *warn = &std::cerr;
#ifdef ITPP_EXCEPTIONS
  static bool it_using_exceptions = true;
#else
  static bool it_using_exceptions = false;
#endif

#endif //DOXYGEN_SHOULD_SKIP_THIS

  void it_assert_f(std::string ass, std::string msg, std::string file, int line)
  {
    char line_str[100];

    std::string error = "*** Assertation failed in ";
    error += file;
    error += " on line ";
    sprintf(line_str, "%d", line);
    error += line_str;
    error += ":\n";
    error += msg;
    error += " (";
    error += ass;
    error += ")";
    std::cerr << error << std::endl << std::flush;
    if (it_using_exceptions)
      throw std::runtime_error(error);
    else
      abort();
  }

  void it_error_f(std::string msg, std::string file, int line)
  {
    char line_str[100];

    std::string error = "*** Error in ";
    error += file;
    error += " on line ";
    sprintf(line_str, "%d", line);
    error += line_str;
    error += ":";
    error += msg;
    std::cerr << error << std::endl << std::flush;
    if (it_using_exceptions)
      throw std::runtime_error(error);
    else
      abort();
  }

  void it_warning_f(std::string msg, std::string file, int line)
  {
    if (warnings_enabled)
      (*warn) << "*** Warning in " << file << " on line " << line << ":" 
	      << std::endl << msg << std::endl << std::flush;
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
