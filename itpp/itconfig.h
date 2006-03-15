/*!
 * \file 
 * \brief Some specific global configurations and definitions
 * \author Tony Ottosson
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

#ifndef ITCONFIG_H
#define ITCONFIG_H

#include <complex>

#ifdef _MSC_VER
#define __WIN32__
#endif

namespace itpp {

  /*! \defgroup parser Argument Parser */
  /*! \defgroup audio Audio */
  /*! \defgroup besselfunctions Bessel Functions */
  /*! \defgroup channels Communication Channel Models */
  /*! \defgroup convertfunc Conversion Functions */
  /*! \defgroup determinant Determinant */
  /*! \defgroup detsource Deterministic Sources */
  /*! \defgroup diag Diagonal Matrices and Functions of Diagonals */
  /*! \defgroup modulators Digital Modulation */
  /*! \defgroup dct Discrete Cosine Transforms (DCT) */
  /*! \defgroup errorfunc Error Functions */
  /*! \defgroup errorhandlingfunc Error Handling Functions */
  /*! \defgroup fft Fast Fourier Transforms (FFT) */
  /*! \defgroup fht Fast Hadamard Transforms */
  /*! \defgroup fastica Fast Independent Component Analysis */
  /*! \defgroup filters Filter Classes and Functions */
  /*! \defgroup fixtypes Fixed-Point Data Types */
  /*! \defgroup fec Forward Error Correcting Codes */
  /*! \defgroup matrix_functions Functions on Matrices */
  /*! \defgroup hypfunc Hyperbolic Functions */
  /*! \defgroup image Image Functions and Classes */
  /*! \defgroup interl Interleavers */
  /*! \defgroup inverse Inverse Matrix */
  /*! \defgroup itfile IT++ File Format */
  /*! \defgroup logexpfunc Logarithmic and Exponential Functions */
  /*! \defgroup lpc LPC-related Functions */
  /*! \defgroup matrixdecomp Matrix Decompositions */
  /*! \defgroup miscfunc Miscellaneous Functions */
  /*! \defgroup integration Numerical Integration */
  /*! \defgroup optimization Numerical Optimization Routines */
  /*! \defgroup poly Polynomial Functions */
  /*! \defgroup randgen Random Number Generation */
  /*! \defgroup reshaping Reshaping of Vectors and Matrices */
  /*! \defgroup sigproc Signal Processing Functions */
  /*! \defgroup linearequations Solving Linear Equation Systems */
  /*! \defgroup sourcecoding Source Coding Routines */
  /*! \defgroup specmat Special Matrices */
  /*! \defgroup statistics Statistics */
  /*! \defgroup timers Timers */
  /*! \defgroup trifunc Trigonometric Functions */
  /*! \defgroup upsample Upsampling of Vectors and Matrices */
  /*! \defgroup windfunc Windowing Functions */

} // namespace itpp


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

#endif // #ifndef ITCONFIG_H
