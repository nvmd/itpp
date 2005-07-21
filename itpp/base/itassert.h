/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2001 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Error handling functions
  \author Tobias Ringström

  $Revision$

  $Date$
*/

#ifndef __itassert_h
#define __itassert_h

#include <itpp/itconfig.h>
#include <iostream>
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

#if ASSERT_LEVEL==0 // No tests
#  define it_assert0(t,s) ((void)0)
#  define it_assert1(t,s) ((void)0)
#elif ASSERT_LEVEL==1 // Only some tests
#  define it_assert0(t,s) ((void)0)
#  define it_assert1(t,s) (void)((t) || (it_assert_f(#t,s,__FILE__,__LINE__),0))
#else // Full tests
  //! Abort if \c t is not true and IT++ is compiled with \c -DASSERT_LEVEL = 2. The message string \c s is printed on standard output.
#  define it_assert0(t,s) (void)((t) || (it_assert_f(#t,s,__FILE__,__LINE__),0))
  //! Abort if \c t is not true and IT++ is compiled with \c -DASSERT_LEVEL = 1 or 2. The message string \c s is printed on standard output.
#  define it_assert1(t,s) (void)((t) || (it_assert_f(#t,s,__FILE__,__LINE__),0))
#endif // ASSERT_LEVEL

  //! Abort if \c t is NOT true and output \c s
#define it_assert(t,s)   (void)((t) || (it_assert_f(#t,s,__FILE__,__LINE__),0))
  //! Abort if \c t is true and output \c s
#define it_error_if(t,s) (void)((!(t)) || (it_error_f(s,__FILE__,__LINE__),0))
  //! Abort and output \c s
#define it_error(s)      it_error_f(s,__FILE__,__LINE__)
  //! Output the warning \c s
#define it_warning(s)    it_warning_f(s,__FILE__,__LINE__)

  //!@}

} //namespace itpp

#endif // __itassert_h
