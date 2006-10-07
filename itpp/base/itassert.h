/*!
 * \file 
 * \brief Definitions of error handling functions
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
 * Copyright (C) 1995-2006  (see AUTHORS file for a list of contributors)
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

#ifndef ITASSERT_H
#define ITASSERT_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <sstream>
#include <string>


namespace itpp {

  /*! 
    \addtogroup errorhandlingfunc

    For all functions, the argument \c s is a string that is displayed.
  
    \code
    it_assert(t,s);   // Abort if \c t is not true.
    it_assert0(t,s);  // Abort if \c t is not true and ASSERT_LEVEL = 2
    it_assert1(t,s);  // Abort if \c t is not true and ASSERT_LEVEL = 2 or 1
    it_error_if(t,s); // Abort if \c t is true.
    it_error(s);      // Abort.
    it_warning(s);    // Show a warning.
    \endcode

    \c it_assert(), \c it_error() and \c it_warning() is always active while
    \c it_assert0() and \c it_assert1() depends on the \c ASSERT_LEVEL variable.
    If \c ASSERT_LEVEL == 0 then none of these are executed while if it is 1
    only \c it_assert1() is executed.
  */
  //!@{
 
  //! Helper function for the \c it_assert functions
  void it_assert_f(std::string ass, std::string msg, std::string file, int line);
  //! Helper function for the \c it_error functions
  void it_error_f(std::string msg, std::string file, int line);
  //! Helper function for the \c it_warning functions
  void it_warning_f(std::string msg, std::string file, int line);

  //! Enable/disable using exceptions for error handling.
  void it_enable_exceptions(bool on);
  //! Enable warnings
  void it_enable_warnings();
  //! Disable warnings
  void it_disable_warnings();
  //! Redirect warnings to the ostream warn_stream
  void it_redirect_warnings(std::ostream *warn_stream);


#define it_assert_base(t,s) \
if(!(t)) { \
  std::ostringstream m_sout; \
  m_sout << s; \
  itpp::it_assert_f(#t,m_sout.str(),__FILE__,__LINE__); \
} else \
  (void)0

#if ASSERT_LEVEL==0 // No tests
#  define it_assert0(t,s) (void)0
#  define it_assert1(t,s) (void)0
#elif ASSERT_LEVEL==1 // Only some tests
#  define it_assert0(t,s) (void)0
#  define it_assert1(t,s) it_assert_base(t,s)
#else // Full tests
  //! Abort if \c t is not true and IT++ is compiled with \c -DASSERT_LEVEL=2
  //! The message string \c s is printed on standard output.
#  define it_assert0(t,s) it_assert_base(t,s)
  //! Abort if \c t is not true and IT++ is compiled with \c -DASSERT_LEVEL=1 
  //! or \c -DASSERT_LEVEL=2. The message string \c s is printed on standard 
  //! output. 
#  define it_assert1(t,s) it_assert_base(t,s)
#endif // ASSERT_LEVEL
  //! Abort if \c t is not true and output \c s
#define it_assert(t,s) it_assert_base(t,s)

  //! Abort if \c t is true and output \c s
#define it_error_if(t,s) \
if((t)) { \
  std::ostringstream m_sout; \
  m_sout << s; \
  itpp::it_error_f(m_sout.str(),__FILE__,__LINE__); \
} else \
  (void)0

  //! Abort and output \c s
#define it_error(s) \
if (true) { \
  std::ostringstream m_sout; \
  m_sout << s; \
  itpp::it_error_f(m_sout.str(),__FILE__,__LINE__); \
} else \
  (void)0

  //! Output the warning \c s
#define it_warning(s) \
if (true) { \
  std::ostringstream m_sout; \
  m_sout << s; \
  itpp::it_warning_f(m_sout.str(),__FILE__,__LINE__); \
} else \
  (void)0

  //!@}

} // namespace itpp

#endif // #ifndef ITASSERT_H
