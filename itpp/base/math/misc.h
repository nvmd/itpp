/*!
 * \file
 * \brief Miscellaneous functions - header file
 * \author Tony Ottosson, Adam Piatyszek and Conrad Sanderson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2008  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#ifndef MISC_H
#define MISC_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <complex>
#include <string>
#include <limits>


namespace std {

  //! \cond

#ifndef HAVE_STD_ISINF
#if (HAVE_DECL_ISINF == 1) || defined(HAVE_ISINF)
  inline int isinf(double x) { return ::isinf(x); }
#elif defined(FPCLASS)
  inline int isinf(double x)
  {
    if (::fpclass(a) == FP_NINF) return -1;
    else if (::fpclass(a) == FP_PINF) return 1;
    else return 0;
  }
#else
  inline int isinf(double x)
  {
    if ((x == x) && ((x - x) != 0.0)) return (x < 0.0 ? -1 : 1);
    else return 0;
  }
#endif // #if (HAVE_DECL_ISINF == 1) || defined(HAVE_ISINF)
#  define HAVE_STD_ISINF 1
#endif // #ifndef HAVE_STD_ISINF

#ifndef HAVE_STD_ISNAN
#if (HAVE_DECL_ISNAN == 1) || defined(HAVE_ISNAN)
  inline int isnan(double x) { return ::isnan(x); }
#else
  inline int isnan(double x) { return ((x != x) ? 1 : 0); }
#endif // #if (HAVE_DECL_ISNAN == 1) || defined(HAVE_ISNAN)
#  define HAVE_STD_ISNAN 1
#endif // #ifndef HAVE_STD_ISNAN

#ifndef HAVE_STD_ISFINITE
#if (HAVE_DECL_ISFINITE == 1) || defined(HAVE_ISFINITE)
  inline int isfinite(double x) { return ::isfinite(x); }
#elif defined(HAVE_FINITE)
  inline int isfinite(double x) { return ::finite(x); }
#else
  inline int isfinite(double x)
  {
    return ((!std::isnan(x) && !std::isinf(x)) ? 1 : 0);
  }
#endif // #if (HAVE_DECL_ISFINITE == 1) || defined(HAVE_ISFINITE)
#  define HAVE_STD_ISFINITE 1
#endif // #ifndef HAVE_STD_ISFINITE

  //! \endcond

  //! Output stream operator for complex numbers
  template <class T>
  std::ostream& operator<<(std::ostream &os, const std::complex<T> &x)
  {
    os << x.real();
    ios::fmtflags saved_format = os.setf(ios::showpos);
    os << x.imag();
    os.setf(saved_format, ios::showpos);
    return os << 'i';
  }

  //! Input stream operator for complex numbers
  template <class T>
  std::istream& operator>>(std::istream &is, std::complex<T> &x)
  {
    T re, im;
    char c;
    is >> c;
    if (c == '(') {
      is >> re >> c;
      if (c == ',') {
        is >> im >> c;
        if (c == ')') {
          x = complex<T>(re, im);
        } else {
          is.setstate(ios_base::failbit);
        }
      } else if (c == ')') {
        x = complex<T>(re, T(0));
      } else {
        is.setstate(ios_base::failbit);
      }
    } else {
      is.putback(c);
      is >> re;
      if (!is.eof() && ((c = is.peek()) == '+' || c == '-')) {
        is >> im >> c;
        if (c == 'i') {
          x = complex<T>(re, im);
        } else {
          is.setstate(ios_base::failbit);
        }
      } else {
        x = complex<T>(re, T(0));
      }
    }
    return is;
  }

} // namespace std


//! itpp namespace
namespace itpp {

  //! Constant Pi
  const double pi = 3.14159265358979323846;

  //! Constant 2*Pi
  const double m_2pi = 2 * pi;

  //! Constant eps
  const double eps = std::numeric_limits<double>::epsilon();

  //! \addtogroup miscfunc
  //!@{

  //! Return true if x is an integer
  inline bool is_int(double x)
  {
    double dummy;
    return (modf(x, &dummy) == 0.0);
  }

  //! Return true if x is an even integer
  inline bool is_even(int x) { return ((x&1) == 0); }


  //! Returns IT++ library version number, e.g. "3.7.1".
  std::string itpp_version();


  //! Returns machine endianness: big-endian = true; little-endian = false
  bool check_big_endianness();

  //!@}

} //namespace itpp


#endif // #ifndef MISC_H
