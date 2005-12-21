/*!
 * \file
 * \brief Definitions of elementary functions on vectors and matrices
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

#ifndef ELMATFUNC_H
#define ELMATFUNC_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/scalfunc.h>
#include <itpp/base/help_functions.h>
#include <itpp/base/converters.h>

namespace itpp {
  // ACTION: This should be split into several files. And documented per group
  /*! \addtogroup logexpfunc */

/*! 
  \brief The gamma function
  \ingroup miscfunc
*/
inline itpp::vec gamma(const itpp::vec &x) {return itpp::vec_function((double(*)(double)) gamma,x);}

/*! 
  \brief The gamma function
  \ingroup miscfunc
*/
inline itpp::mat gamma(const itpp::mat &x) {return itpp::mat_function((double(*)(double)) gamma,x);}

  //!\addtogroup logexpfunc
  //!@{

  //! Exp of the elements
  inline vec exp(const vec &x) {return vec_function((double(*)(double)) std::exp,x);}
  //! Exp of cvec
  inline cvec exp(const cvec &x) {return cvec_function((std::complex<double>(*)(const std::complex<double> &)) std::exp,x);}

  //! Exp of the elements
  inline mat exp(const mat &x) {return mat_function((double(*)(double)) std::exp,x);}
  //! Exp of cmat
  inline cmat exp(const cmat &x) {return cmat_function((std::complex<double>(*)(const std::complex<double> &)) std::exp,x);}
  //! Calculates x to the power of y (x^y)
  inline vec pow(const double x, const vec &y)
    {return double_vec_function((double(*)(double,double)) std::pow,x,y);}
  //! Calculates x to the power of y (x^y)
  inline mat pow(const double x, const mat &y)
    {return double_mat_function((double(*)(double,double)) std::pow,x,y);}
  //! Calculates x to the power of y (x^y)
  inline vec pow(const vec &x, const double y)
    {return vec_double_function((double(*)(double,double)) std::pow,x,y);}
  //! Calculates x to the power of y (x^y)
  inline mat pow(const mat &x, const double y)
    {return mat_double_function((double(*)(double,double)) std::pow,x,y);}
  //! Calculates two to the power of x (2^x)
  inline vec pow2(const vec &x)
    {return vec_function((double(*)(double)) pow2,x);}
  //! Calculates two to the power of x (2^x)
  inline mat pow2(const mat &x)
    {return mat_function((double(*)(double)) pow2,x);}
  //! Calculates ten to the power of x (10^x)
  inline vec pow10(const vec &x)
    {return vec_function((double(*)(double)) pow10,x);}
  //! Calculates ten to the power of x (10^x)
  inline mat pow10(const mat &x)
    {return mat_function((double(*)(double)) pow10,x);}

  //! The natural logarithm of the elements
  inline vec log(const vec &x) {return vec_function((double(*)(double)) std::log,x);}
  //! The natural logarithm of the elements
  inline mat log(const mat &x) {return mat_function((double(*)(double)) std::log,x);}

  //! log-2 of the elements
  inline vec log2(const vec &x) {return vec_function((double(*)(double)) itpp::log2,x);}
  //! log-2 of the elements
  inline mat log2(const mat &x) {return mat_function((double(*)(double)) itpp::log2,x);}

  //! log-10 of the elements
  inline vec log10(const vec &x) {return vec_function((double(*)(double)) std::log10,x);}
  //! log-10 of the elements
  inline mat log10(const mat &x) {return mat_function((double(*)(double)) std::log10,x);}
  //! log-b of \c x
  inline vec logb(const short b,const vec &x) {return double_vec_function(itpp::logb,b,x);}
  //! log-b of \c x
  inline mat logb(const short b,const mat &x) {return double_mat_function(itpp::logb,b,x);}

  //! Calculates 10*log10(x)
  inline vec dB(const vec &x) {return vec_function(dB,x);}
  //! Calculates 10*log10(x)
  inline mat dB(const mat &x) {return mat_function(dB,x);}
  //! Calulates the inverse of dB, 10^(x/10)
  inline vec inv_dB(const vec &x) {return vec_function(inv_dB,x);}
  //! Calculates the inverse of dB, 10^(x/10)
  inline mat inv_dB(const mat &x) {return mat_function(inv_dB,x);}
  //!@}


  //!\addtogroup errorfunc
  //!@{

  //! Error function
  inline vec erf(const vec &x) {return vec_function((double(*)(double)) ::erf,x);}
  //! Error function
  inline mat erf(const mat &x) {return mat_function((double(*)(double)) ::erf,x);}
  //! Error function
  inline cvec erf(const cvec &x) {return cvec_function((std::complex<double>(*)(const std::complex<double> &)) erf,x);}
  //! Error function
  inline cmat erf(const cmat &x) {return cmat_function((std::complex<double>(*)(const std::complex<double> &)) erf,x);}
  //! Inverse of error function
  inline vec erfinv(const vec &x) {return vec_function(erfinv,x);}
  //! Inverse of error function
  inline mat erfinv(const mat &x) {return mat_function(erfinv,x);}
  //! Complementary error function
  inline vec erfc(const vec &x) {return vec_function((double(*)(double)) ::erfc,x);}
  //!Complementary error function
  inline mat erfc(const mat &x) {return mat_function((double(*)(double)) ::erfc,x);}
  //! Q-function
  inline vec Qfunc(const vec &x) {return vec_function(Qfunc,x);}
  //! Q-function
  inline mat Qfunc(const mat &x) {return mat_function(Qfunc,x);}
  //!@}


  //!\addtogroup trifunc
  //!@{

  //! Sine function
  inline vec sin(const vec &x) {return vec_function((double(*)(double)) std::sin,x);}
  //! Sine function
  inline mat sin(const mat &x) {return mat_function((double(*)(double)) std::sin,x);}
  //! Cosine function
  inline vec cos(const vec &x) {return vec_function((double(*)(double)) std::cos,x);}
  //! Cosine function
  inline mat cos(const mat &x) {return mat_function((double(*)(double)) std::cos,x);}
  //! Tan function
  inline vec tan(const vec &x) {return vec_function((double(*)(double)) std::tan,x);}
  //! Tan function
  inline mat tan(const mat &x) {return mat_function((double(*)(double)) std::tan,x);}
  //! Inverse sine function
  inline vec asin(const vec &x) {return vec_function((double(*)(double)) std::asin,x);}
  //! Inverse sine function
  inline mat asin(const mat &x) {return mat_function((double(*)(double)) std::asin,x);}
  //! Inverse cosine function
  inline vec acos(const vec &x) {return vec_function((double(*)(double)) std::acos,x);}
  //! Inverse cosine function
  inline mat acos(const mat &x) {return mat_function((double(*)(double)) std::acos,x);}
  //! Inverse tan function
  inline vec atan(const vec &x) {return vec_function((double(*)(double)) std::atan,x);}
  //! Inverse tan function
  inline mat atan(const mat &x) {return mat_function((double(*)(double)) std::atan,x);}
  //! Sinc function, sin(pi*x)/(pi*x)
  inline vec sinc(const vec &x) {return vec_function((double(*)(double)) sinc,x);}
  //! Sinc function, sin(pi*x)/(pi*x)
  inline mat sinc(const mat &x) {return mat_function((double(*)(double)) sinc,x);}
  //!@}


  //!\addtogroup hypfunc
  //!@{

  //! Sine hyperbolic function
  inline vec sinh(const vec &x) {return vec_function((double(*)(double)) std::sinh,x);}
  //! Sine hyperbolic function
  inline mat sinh(const mat &x) {return mat_function((double(*)(double)) std::sinh,x);}
  //! Cosine hyperbolic function
  inline vec cosh(const vec &x) {return vec_function((double(*)(double)) std::cosh,x);}
  //! Cosine hyperbolic function
  inline mat cosh(const mat &x) {return mat_function((double(*)(double)) std::cosh,x);}
  //! Tan hyperbolic function
  inline vec tanh(const vec &x) {return vec_function((double(*)(double)) std::tanh,x);}
  //! Tan hyperbolic function
  inline mat tanh(const mat &x) {return mat_function((double(*)(double)) std::tanh,x);}
  //! Inverse sine hyperbolic function
  inline vec asinh(const vec &x) {return vec_function((double(*)(double)) ::asinh,x);}
  //! Inverse sine hyperbolic function
  inline mat asinh(const mat &x) {return mat_function((double(*)(double)) ::asinh,x);}
  //! Inverse cosine hyperbolic function
  inline vec acosh(const vec &x) {return vec_function((double(*)(double)) ::acosh,x);}
  //! Inverse cosine hyperbolic function
  inline mat acosh(const mat &x) {return mat_function((double(*)(double)) ::acosh,x);}
  //! Inverse tan hyperbolic function
  inline vec atanh(const vec &x) {return vec_function((double(*)(double)) ::atanh,x);}
  //! Inverse tan hyperbolic function
  inline mat atanh(const mat &x) {return mat_function((double(*)(double)) ::atanh,x);}
  //!@}


  //!\addtogroup miscfunc
  //!@{

  //! Round to nearest upper integer
  inline vec ceil(const vec &x) {return vec_function((double(*)(double)) std::ceil,x);}
  //! Round to nearest upper integer
  inline mat ceil(const mat &x) {return mat_function((double(*)(double)) std::ceil,x);}
  //! Round to nearest lower integer
  inline vec floor(const vec &x) {return vec_function((double(*)(double)) std::floor,x);}
  //! Round to nearest lower integer
  inline mat floor(const mat &x) {return mat_function((double(*)(double)) std::floor,x);}

  //! Round to nearest integer
  inline vec round(const vec &x) {return vec_function((double(*)(double)) round,x);}
  //! Round to nearest integer
  inline mat round(const mat &x) {return mat_function((double(*)(double)) round,x);}

  //! Round to nearest integer and return ivec
  inline ivec round_i(const vec &x) {return to_ivec(round(x));}
  //! Round to nearest integer and return imat
  inline imat round_i(const mat &x) {return to_imat(round(x));}


  //! Absolute value
  inline vec abs(const vec &x) {return vec_function((double(*)(double)) std::fabs,x);}
  //! Absolute value
  inline mat abs(const mat &x) {return mat_function((double(*)(double)) std::fabs,x);}
  //! Absolute value
  ivec abs(const ivec &x);
  //! Absolute value
  imat abs(const imat &x);
  //! Square of elements
  inline vec sqr(const vec &x) {return vec_function(sqr,x);}
  //! Square of elements
  inline mat sqr(const mat &x) {return mat_function(sqr,x);}
  //! Square of elements
  vec sqr(const cvec &x);
  //! Square of elements
  mat sqr(const cmat &x);
  //! Signum function
  inline vec sign(const vec &x) {return vec_function(sign,x);}
  //! Signum function
  inline mat sign(const mat &x) {return mat_function(sign,x);}

  //! Square root of the elements
  inline vec sqrt(const vec &x) {return vec_function((double(*)(double)) std::sqrt,x);}
  //! Square root of the elements
  inline mat sqrt(const mat &x) {return mat_function((double(*)(double)) std::sqrt,x);}

  // ACTION: mod()

  //! Elementwise reminder of the division x/y for vec and double
  inline vec rem(const vec &x,const double &y) {return vec_double_function(rem,x,y);}
  //! Elementwise reminder of the division x/y for double and vec
  inline vec rem(const double &x,const vec &y) {return double_vec_function(rem,x,y);}
  //! Elementwise reminder of the division x/y for mat and double
  inline mat rem(const mat &x,const double &y) {return mat_double_function(rem,x,y);}
  //! Elementwise reminder of the division x/y for double and mat
  inline mat rem(const double &x,const mat &y) {return double_mat_function(rem,x,y);}

  //! Absolute value
  vec abs(const cvec &x);
  //! Absolute value
  mat abs(const cmat &x);
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
  cvec conj(const cvec &x);
  //! Conjugate of complex value
  cmat conj(const cmat &x);

  //! Returns \a true if all elements are ones and \a false otherwise
  bool all(const Vec<bin> &testvec);
  //! Returns \a true if any element is one and \a false otherwise
  bool any(const Vec<bin> &testvec);

  //! Round each element to zero if element < threshold
  inline vec round_to_zero(const vec &x, double threshold = 1e-15) {
    return vec_double_function((double(*)(double,double)) round_to_zero,x,threshold);
  }
  //! Round each element to zero if element < threshold
  inline mat round_to_zero(const mat &x, double threshold = 1e-15) {
    return mat_double_function((double(*)(double,double)) round_to_zero,x,threshold);
  }
  //! Round each element to zero if element < threshold
  cvec round_to_zero(const cvec &x, double threshold = 1e-15);

  //! Round each element to zero if element < threshold
  cmat round_to_zero(const cmat &x, double threshold = 1e-15);

  //!@}

} // namespace itpp

#endif // #ifndef ELMATFUNC_H




