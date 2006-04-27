/*!
 * \file
 * \brief Elementary mathematical functions
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

#ifndef ELMATFUNC_H
#define ELMATFUNC_H

#include <itpp/base/help_functions.h>
#include <itpp/base/converters.h>
#include <cmath>


#ifdef _MSC_VER

// These functions are part of C99. But apparently not in Visual C++.

//!\addtogroup miscfunc
//!@{

//! True gamma function
double tgamma(double x);

//! Lograrithm of an absolute gamma function
extern int signgam;
double lgamma(double x);

//! Cubic root
double cbrt(double x);
//!@}

#endif

namespace itpp {

  //!\addtogroup miscfunc
  //!@{

  // -------------------- sqr function --------------------

  //! Square of x
  inline double sqr(double x) { return (x * x); }
  //! Square of complex-valued x, ||x||^2
  inline double sqr(const std::complex<double>& x) 
  { 
    return (x.real() * x.real() + x.imag() * x.imag()); 
  }
  //! Square of elements
  inline vec sqr(const vec &x) { return apply_function<double>(sqr, x); }
  //! Square of elements
  inline mat sqr(const mat &x) { return apply_function<double>(sqr, x); }
  //! Square of elements
  vec sqr(const cvec &x); 
  //! Square of elements
  mat sqr(const cmat &x);


  // -------------------- abs function --------------------

//  //! Absolute value
//  inline signed char abs(signed char x) { return (x >= 0 ? x : -x); }
  //! Absolute value
  inline short abs(short x) { return (x >= 0 ? x : -x); }
  //! Absolute value
  inline int abs(int x) { return (x >= 0 ? x : -x); }
  //! Absolute value
  inline vec abs(const vec &x) { return apply_function<double>(std::fabs, x); }
  //! Absolute value
  inline mat abs(const mat &x) { return apply_function<double>(std::fabs, x); }
  //! Absolute value
  inline ivec abs(const ivec &x) { return apply_function<int>(abs, x); }
  //! Absolute value
  inline imat abs(const imat &x) { return apply_function<int>(abs, x); }
  //! Absolute value
  vec abs(const cvec &x);
  //! Absolute value
  mat abs(const cmat &x);


  // -------------------- sign/sgn functions --------------------
  
  //! The sign of x
  inline double sign(double x)
  { 
    return (x == 0.0 ? 0.0 : (x < 0.0 ? -1.0 : 1.0));
  }
  //! Signum function
  inline vec sign(const vec &x) { return apply_function<double>(sign, x); }
  //! Signum function
  inline mat sign(const mat &x) { return apply_function<double>(sign, x); }
  //! The sign of x
  inline double sgn(double x)
  { 
    return (x == 0.0 ? 0.0 : (x < 0.0 ? -1.0 : 1.0));
  }
  //! Signum function
  inline vec sgn(const vec &x) { return apply_function<double>(sgn, x); }
  //! Signum function
  inline mat sgn(const mat &x) { return apply_function<double>(sgn, x); }


  // -------------------- sqrt function --------------------

  //! Square root of the elements
  inline vec sqrt(const vec &x) { return apply_function<double>(std::sqrt, x); }
  //! Square root of the elements
  inline mat sqrt(const mat &x) { return apply_function<double>(std::sqrt, x); }


  // -------------------- gamma function --------------------

  //! Gamma function
  inline double gamma(double x) 
  { 
    // lgamma() needs to be executed before, since it sets signgam
    double lg = lgamma(x);
    return signgam * std::exp(lg);
  }
  //! The gamma function
  inline vec gamma(const vec &x) { return apply_function<double>(gamma, x); }
  //! The gamma function
  inline mat gamma(const mat &x) { return apply_function<double>(gamma, x); }

  
  // -------------------- rem function --------------------

  //! The reminder of the division x/y
  inline double rem(double x, double y) { return fmod(x, y); }
  //! Elementwise reminder of the division x/y for vec and double
  inline vec rem(const vec &x, double y) 
  { 
    return apply_function<double>(rem, x, y); 
  }
  //! Elementwise reminder of the division x/y for double and vec
  inline vec rem(double x, const vec &y) 
  {
    return apply_function<double>(rem, x, y);
  }
  //! Elementwise reminder of the division x/y for mat and double
  inline mat rem(const mat &x, double y) 
  {
    return apply_function<double>(rem, x, y);
  }
  //! Elementwise reminder of the division x/y for double and mat
  inline mat rem(double x, const mat &y) 
  {
    return apply_function<double>(rem, x, y);
  }

  // -------------------- mod function --------------------

  //! Calculates the modulus, i.e. the signed reminder after division
  inline int mod(int k, int n)
  {
    return (n == 0) ? k : (k - n * floor_i(static_cast<double>(k) / n ));
  }


  // -------------------- complex related functions --------------------

  //! Real part of complex values
  vec real(const cvec &x);
  //! Real part of complex values
  mat real(const cmat &x);
  //! Imaginary part of complex values
  vec imag(const cvec &x);
  //! Imaginary part of complex values
  mat imag(const cmat &x);

  //! Argument (angle)
  vec arg(const cvec &x);
  //! Argument (angle)
  mat arg(const cmat &x);
  //! Angle
  inline vec angle(const cvec &x) { return arg(x); }
  //! Angle
  inline mat angle(const cmat &x) { return arg(x); }

  //! Conjugate of complex value
  inline cvec conj(const cvec &x)
  { 
    return apply_function<std::complex<double> >(std::conj, x); 
  }
  //! Conjugate of complex value
  inline cmat conj(const cmat &x)
  { 
    return apply_function<std::complex<double> >(std::conj, x); 
  }

  //!@}

} // namespace itpp

#endif // #ifndef ELMATFUNC_H




