/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2004 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Some it++ specific configurations and definitions

  $Revision$

  $Date$
*/

#ifndef __itconfig_h
#define __itconfig_h

#include <complex>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define ITPP_DEFAULT_EXCEPTIONS 0
#endif //DOXYGEN_SHOULD_SKIP_THIS


#ifdef _MSC_VER
#define __WIN32__
#endif

//---------------------------------------------------------------

//! Output stream operator for complex numbers
namespace std {
  template <class T>
    std::ostream& operator<<(std::ostream &os, const std::complex<T> &x)
    {
      os <<  x.real() ;
      if (x.imag() >= 0) os << '+' << x.imag();
      else os << x.imag();
      return os << 'i';
    }
}

//! Input stream operator for complex numbers
namespace std {
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
}

#endif // __itconfig_h
