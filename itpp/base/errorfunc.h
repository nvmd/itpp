/*!
 * \file
 * \brief Error functions - header file
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

#ifndef ERRORFUNC_H
#define ERRORFUNC_H

#include <itpp/base/help_functions.h>
//#include <cmath>


#ifdef _MSC_VER

// These functions are part of C99. But apparently not in Visual C++.

//!\addtogroup errorfunc
//!@{

//! \brief Error function
double erf(double x);

//! Complementary error function
double erfc(double x);
//!@}

#endif

namespace itpp {

  //!\addtogroup errorfunc
  //!@{

  /*!
   * \brief Error function for complex argument
   * \author Adam Piatyszek
   *
   * This function calculates a well known error function \c erf(z)
   * for complex \c z. The implementation is based on unofficial
   * implementation for Octave. Here is a part of the author's note
   * from original sources:
   *
   * Put together by John Smith john at arrows dot demon dot co dot uk, 
   * using ideas by others.
   *
   * Calculate \c erf(z) for complex \c z.
   * Three methods are implemented; which one is used depends on z.
   *
   * The code includes some hard coded constants that are intended to
   * give about 14 decimal places of accuracy. This is appropriate for
   * 64-bit floating point numbers. 
   */
  std::complex<double> erf(const std::complex<double>& z);

  //! Inverse of error function
  double erfinv(double x);

  //! Q-function
  double Qfunc(double x);


  // ----------------------------------------------------------------------
  // functions for matrices and vectors
  // ----------------------------------------------------------------------

  //! Error function
  inline vec erf(const vec &x) 
  {
    return apply_function<double>(::erf, x);
  }
  //! Error function
  inline mat erf(const mat &x) 
  {
    return apply_function<double>(::erf, x);
  }
  //! Error function
  inline cvec erf(const cvec &x) 
  {
    return apply_function<std::complex<double> >(erf, x);
  }
  //! Error function
  inline cmat erf(const cmat &x) 
  { 
    return apply_function<std::complex<double> >(erf, x);
  }

  //! Inverse of error function
  inline vec erfinv(const vec &x) 
  {
    return apply_function<double>(erfinv, x);
  }
  //! Inverse of error function
  inline mat erfinv(const mat &x) 
  {
    return apply_function<double>(erfinv, x);
  }

  //! Complementary error function
  inline vec erfc(const vec &x) 
  {
    return apply_function<double>(::erfc, x);
  }
  //! Complementary error function
  inline mat erfc(const mat &x) 
  {
    return apply_function<double>(::erfc, x);
  }

  //! Q-function
  inline vec Qfunc(const vec &x) 
  {
    return apply_function<double>(Qfunc, x);
  }

  //! Q-function
  inline mat Qfunc(const mat &x) 
  {
    return apply_function<double>(Qfunc, x);
  }
  //!@}

} // namespace itpp

#endif // #ifndef ERRORFUNC_H




