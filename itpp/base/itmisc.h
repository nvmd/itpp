/*!
 * \file
 * \brief Definition of IT++ miscellaneous functions
 * \author Tony Ottosson, Adam Piatyszek and Conrad Sanderson
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

#ifndef ITMISC_H
#define ITMISC_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <complex>
#include <string>
#include <limits>


namespace std {

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


//!\addtogroup miscfunc
//!@{

#ifndef HAVE_ISNAN
/*!
 * \brief Check if \c x is NaN (Not a Number)
 * \note Emulation of a C99 function via the IEEE 754 standard 
 */
inline int isnan(double x)
{ 
  if (x != x) return 1;
  else return 0;
}
#endif

#ifndef HAVE_ISINF
/*!
 * \brief Check if \c x is either -Inf or +Inf
 *
 * Returns -1 if \c x is -Inf or +1 if \c x is +Inf. Otherwise returns 0.
 * 
 * \note Emulation of a C99 function via the IEEE 754 standard 
 */
inline int isinf(double x) 
{ 
  if ((x == x) && ((x - x) != 0.0)) 
    return (x < 0.0 ? -1 : 1);
  else return 0;
}
#endif

#ifndef HAVE_FINITE
/*!
 * \brief Check if \c x is a finite floating point number
 * \note Emulation of a C99 function via the IEEE 754 standard 
 */
inline int finite(double x) 
{ 
  if (!isnan(x) && !isinf(x)) return 1;
  else return 0;
}
#endif

#ifndef HAVE_ISFINITE
/*!
 * \brief Check if \c x is a finite floating point number
 * \note Emulation of a C99 function via the IEEE 754 standard 
 */
inline int isfinite(double x) { return finite(x); }
#endif

//!@}


#endif // #ifndef ITMISC_H
