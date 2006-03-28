/*!
 * \file
 * \brief Definition of IT++ miscleaneous functions
 * \author Tony Ottosson and Adam Piatyszek
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

#ifndef ITMISC_H
#define ITMISC_H

#include <complex>
#include <string>


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

  /*!
   * \brief Return IT++ version
   *
   * Returns the version number of the IT++ library, e.g. "3.7.1".
   */
  std::string itpp_version(void);

} //namespace itpp


#endif // #ifndef ITMISC_H
