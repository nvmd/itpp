/*!
 * \file
 * \brief Elementary mathematical functions - header file
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

#ifndef ELEM_MATH_H
#define ELEM_MATH_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <itpp/base/help_functions.h>
#include <itpp/base/converters.h>
#include <cstdlib>


//!\addtogroup miscfunc
//!@{

#ifndef HAVE_TGAMMA
//! True gamma function
double tgamma(double x);
#endif

#if !defined(HAVE_LGAMMA) || (HAVE_DECL_SIGNGAM != 1)
//! Lograrithm of an absolute gamma function
double lgamma(double x);
//! Global variable needed by \c lgamma function
extern int signgam;
#endif

#ifndef HAVE_CBRT
//! Cubic root
double cbrt(double x);
#endif

//!@}

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

  //! Absolute value
  inline vec abs(const vec &x) { return apply_function<double>(std::fabs, x); }
  //! Absolute value
  inline mat abs(const mat &x) { return apply_function<double>(std::fabs, x); }
  //! Absolute value
  inline ivec abs(const ivec &x) { return apply_function<int>(std::abs, x); }
  //! Absolute value
  inline imat abs(const imat &x) { return apply_function<int>(std::abs, x); }
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
  inline double gamma(double x) { return tgamma(x); }
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


  // -------------------- factorial coefficient function --------------------

  //! Calculates factorial coefficient for index <= 170.
  double fact(int index);


  // -------------------- binomial coefficient function --------------------

  //! Compute the binomial coefficient "n over k".
  double binom(int n, int k);

  //! Compute the binomial coefficient "n over k".
  int binom_i(int n, int k);

  //! Compute the base 10 logarithm of the binomial coefficient "n over k".
  double log_binom(int n, int k);


  // -------------------- greatest common divisor function --------------------

  /*!
   * \brief Compute the greatest common divisor (GCD) \a g of the elements
   * \a a and \a b.
   *
   * \a a and \a b must be non-negative integers. \a gdc(0, 0) is 0 by
   * convention; all other GCDs are positive integers.
   */
  int gcd(int a, int b);


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

  // Added due to a failure in MSVC++ .NET 2005, which crashes on this
  // code.
#ifndef _MSC_VER
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
#else
  //! Conjugate of complex value
  cvec conj(const cvec &x);

  //! Conjugate of complex value
  cmat conj(const cmat &x);
#endif

  //!@}

} // namespace itpp

#endif // #ifndef ELEM_MATH_H




