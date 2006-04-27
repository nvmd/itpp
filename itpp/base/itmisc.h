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

  //! Constant pi
  const double pi = 3.14159265358979323846;

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


  //! Compute the binomial coefficient "n over k" as a float.
  double binom(int n, int k);

  //! Compute the binomial coefficient "n over k" as an integer.
  int binom_i(int n, int k);

  //! Compute the base 10 logarithm of the binomial coefficient "n over k".
  double log_binom(int n, int k);


  //double sigmoid(double x) { return 1.0/(1.0+exp(-x)); }

  //! Calculates factorial coefficient for index <= 170.
  double fact(int index);

  /*!
    \brief Returns the greatest common divisor (GCD) \a g of the elements \a
           a and \a b. 

    \a a and \a b must be non-negative integers. \a gdc(0,0) is 0 by
    convention; all other GCDs are positive integers.
  */
  long gcd(long a, long b);


  /*!
   * \brief Return IT++ version
   *
   * Returns the version number of the IT++ library, e.g. "3.7.1".
   */
  std::string itpp_version(void);

  //!@}

} //namespace itpp


#endif // #ifndef ITMISC_H
