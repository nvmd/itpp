/*!
 * \file
 * \brief Logarithmic and exponenential functions - header file
 * \author Tony Ottosson, Adam Piatyszek and Conrad Sanderson
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

#ifndef LOG_EXP_H
#define LOG_EXP_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <itpp/base/help_functions.h>
#include <limits>


/*!
 * \addtogroup logexpfunc
 * @{
 */

#ifndef HAVE_LOG1P
//! Lograrithm of an argument \c x plus one
inline double log1p(double x) { return std::log(1.0 + x); }
#endif

#ifndef HAVE_LOG2
#undef log2                     // This is required at least for Cygwin
//! Base-2 logarithm
inline double log2(double x)
{
  return (std::log(x) * 1.442695040888963387004650940070860087871551513671875);
}
#endif

/*!
 * @}
 */


namespace itpp {

  //!\addtogroup logexpfunc
  //!@{

  // ----------------------------------------------------------------------
  // scalar functions
  // ----------------------------------------------------------------------

  //! Base-b logarithm
  inline double logb(double b, double x)
  {
    return (std::log(x) / std::log(b));
  }

  //! Calculate two to the power of x (2^x); x is integer
  inline int pow2i(int x) { return ((x < 0) ? 0 : (1 << x)); }
  //! Calculate two to the power of x (2^x)
  inline double pow2(double x) { return pow(2.0, x); }

  //! Calculate ten to the power of x (10^x)
  inline double pow10(double x) { return pow(10.0, x); }

  //! Decibel of x (10*log10(x))
  inline double dB(double x) { return 10.0 * log10(x); }
  //! Inverse of decibel of x
  inline double inv_dB(double x) { return pow(10.0, 0.1 * x); }

  //! Calculate the number of bits needed to represent an inteager \c n
  inline int int2bits(int n)
  {
    it_assert(n >= 0, "int2bits(): Improper argument value");

    if (n == 0)
      return 1;

    int b = 0;
    while (n) {
      n >>= 1;
      ++b;
    }
    return b;
  }

  //! Calculate the number of bits needed to represent \c n different values (levels).
  inline int levels2bits(int n)
  {
    it_assert(n > 0,"levels2bits(): Improper argument value");
    return int2bits(--n);
  }

  //! Deprecated function. Please use int2bits() or levels2bits() instead.
  inline int needed_bits(int n)
  {
    it_warning("needed_bits(): This function is depreceted. Depending on your needs, please use int2bits() or levels2bits() instead.");
    return int2bits(n);
  }

  //! Constant definition to speed up trunc_log() and trunc_exp()
  const double log_double_max = std::log(std::numeric_limits<double>::max());
  //! Constant definition to speed up trunc_log(), trunc_exp() and log_add()
  const double log_double_min = std::log(std::numeric_limits<double>::min());

  /*!
    \brief Truncated natural logarithm function

    This truncated function provides a solution in the cases when the
    logarithm argument is less or equal to zero or infinity. The function
    checks for such extreme values and use some kind of truncation
    (saturation) before calculating the logarithm.

    The truncated logarithm function can be used for calculation of
    log-likelihood in soft demodulators, when numerical instability problem
    might occur.
  */
  inline double trunc_log(double x)
  {
    if (std::numeric_limits<double>::is_iec559) {
      if (x == std::numeric_limits<double>::infinity())
	return log_double_max;
      if (x <= 0)
	return log_double_min;
    }
    return std::log(x);
  }

  /*!
    \brief Truncated exponential function

    This truncated function provides a solution in the case when the
    exponent function results in infinity. The function checks for an
    extreme value and use truncation (saturation) before calculating
    the result.

    The truncated exponential function can be used  when numerical
    instability problem occurs for a standard exp function.
  */
  inline double trunc_exp(double x)
  {
    if (std::numeric_limits<double>::is_iec559
	&& (x >= log_double_max))
      return std::numeric_limits<double>::max();
    return std::exp(x);
  }


  //! Safe substitute for <tt>log(exp(log_a) + exp(log_b))</tt>
  inline double log_add(double log_a, double log_b)
  {
    if (log_a < log_b) {
      double tmp = log_a;
      log_a = log_b;
      log_b = tmp;
    }
    double negdelta = log_b - log_a;
    if (negdelta < log_double_min)
      return log_a;
    else
      return (log_a + log1p(std::exp(negdelta)));
  }


  // ----------------------------------------------------------------------
  // functions on vectors and matrices
  // ----------------------------------------------------------------------

  //! Exp of the elements of a vector \c x
  inline vec exp(const vec &x)
  {
    return apply_function<double>(std::exp, x);
  }
  //! Exp of the elements of a complex vector \c x
  inline cvec exp(const cvec &x)
  {
    return apply_function<std::complex<double> >(std::exp, x);
  }
  //! Exp of the elements of a matrix \c m
  inline mat exp(const mat &m)
  {
    return apply_function<double>(std::exp, m);
  }
  //! Exp of the elements of a complex matrix \c m
  inline cmat exp(const cmat &m)
  {
    return apply_function<std::complex<double> >(std::exp, m);
  }

  //! Calculates x to the power of y (x^y)
  inline vec pow(const double x, const vec &y)
  {
    return apply_function<double>(std::pow, x, y);
  }
  //! Calculates x to the power of y (x^y)
  inline mat pow(const double x, const mat &y)
  {
    return apply_function<double>(std::pow, x, y);
  }
  //! Calculates x to the power of y (x^y)
  inline vec pow(const vec &x, const double y)
  {
    return apply_function<double>(std::pow, x, y);
  }
  //! Calculates x to the power of y (x^y)
  inline mat pow(const mat &x, const double y)
  {
    return apply_function<double>(std::pow, x, y);
  }

  //! Calculates two to the power of x (2^x)
  inline vec pow2(const vec &x)
  {
    return apply_function<double>(pow2, x);
  }
  //! Calculates two to the power of x (2^x)
  inline mat pow2(const mat &x)
  {
    return apply_function<double>(pow2, x);
  }

  //! Calculates ten to the power of x (10^x)
  inline vec pow10(const vec &x)
  {
    return apply_function<double>(pow10, x);
  }
  //! Calculates ten to the power of x (10^x)
  inline mat pow10(const mat &x)
  {
    return apply_function<double>(pow10, x);
  }

  //! The natural logarithm of the elements
  inline vec log(const vec &x)
  {
    return apply_function<double>(std::log, x);
  }
  //! The natural logarithm of the elements
  inline mat log(const mat &x)
  {
    return apply_function<double>(std::log, x);
  }
  //! The natural logarithm of the elements
  inline cvec log(const cvec &x)
  {
    return apply_function<std::complex<double> >(std::log, x);
  }
  //! The natural logarithm of the elements
  inline cmat log(const cmat &x)
  {
    return apply_function<std::complex<double> >(std::log, x);
  }

  //! log-2 of the elements
  inline vec log2(const vec &x)
  {
    return apply_function<double>(::log2, x);
  }
  //! log-2 of the elements
  inline mat log2(const mat &x)
  {
    return apply_function<double>(::log2, x);
  }

  //! log-10 of the elements
  inline vec log10(const vec &x)
  {
    return apply_function<double>(std::log10, x);
  }
  //! log-10 of the elements
  inline mat log10(const mat &x)
  {
    return apply_function<double>(std::log10, x);
  }

  //! log-b of \c x
  inline vec logb(double b, const vec &x)
  {
    return apply_function<double>(itpp::logb, b, x);
  }
  //! log-b of \c x
  inline mat logb(double b, const mat &x)
  {
    return apply_function<double>(itpp::logb, b, x);
  }

  //! Calculates 10*log10(x)
  inline vec dB(const vec &x)
  {
    return apply_function<double>(dB, x);
  }
  //! Calculates 10*log10(x)
  inline mat dB(const mat &x)
  {
    return apply_function<double>(dB, x);
  }

  //! Calulates the inverse of dB, 10^(x/10)
  inline vec inv_dB(const vec &x)
  {
    return apply_function<double>(inv_dB, x);
  }
  //! Calculates the inverse of dB, 10^(x/10)
  inline mat inv_dB(const mat &x)
  {
    return apply_function<double>(inv_dB, x);
  }

  //! Deprecated function. Please use int2bits() or levels2bits() instead.
  inline ivec needed_bits(const ivec& v)
  {
    it_warning("needed_bits(): This function is depreceted. Depending on your needs, please use int2bits() or levels2bits() instead.");
    return apply_function<int>(int2bits, v);
  }

  //! Calculate the number of bits needed to represent each inteager in a vector
  inline ivec int2bits(const ivec& v)
  {
    return apply_function<int>(int2bits, v);
  }

  //! Calculate the number of bits needed to represent a numer of levels saved in a vector
  inline ivec levels2bits(const ivec& v)
  {
    return apply_function<int>(levels2bits, v);
  }

  //!@}

} // namespace itpp

#endif // #ifndef LOG_EXP_H




