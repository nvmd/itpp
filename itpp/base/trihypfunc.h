/*!
 * \file
 * \brief Trigonometric and hyperbolic functions
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

#ifndef TRIHYPFUNC_H
#define TRIHYPFUNC_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <itpp/base/help_functions.h>


//! \addtogroup hypfunc
//!@{

#ifndef HAVE_ASINH
//! Arcus sinhyp
inline double asinh(double x) 
{
  return ((x >= 0) ? log(x + sqrt(x * x + 1)) : -log( -x + sqrt(x * x + 1)));
}
#endif

#ifndef HAVE_ACOSH
//! Arcus coshyp
inline double acosh(double x)
{
  it_error_if(x < 1, "acosh(): Argument must be greater then 1.");
  return log(x + sqrt(x * x - 1.0));
}
#endif

#ifndef HAVE_ATANH
//! Arcus tanhyp
inline double atanh(double x)
{
  it_error_if(fabs(x) >= 1,"atanh(): Argument out of range.");
  return 0.5 * log((x + 1) / (x - 1));
}
#endif

//!@}


namespace itpp {

  //!\addtogroup trifunc
  //!@{

  //! Sinc function: sinc(x) = sin(pi*x)/pi*x
  inline double sinc(double x)
  { 
    return ((x == 0) ? 1.0 : sin(itpp::pi * x) / itpp::pi / x); 
  }

  //! Sine function
  inline vec sin(const vec &x) { return apply_function<double>(std::sin, x); }
  //! Sine function
  inline mat sin(const mat &x) { return apply_function<double>(std::sin, x); }
  //! Cosine function
  inline vec cos(const vec &x) { return apply_function<double>(std::cos, x); }
  //! Cosine function
  inline mat cos(const mat &x) { return apply_function<double>(std::cos, x); }
  //! Tan function
  inline vec tan(const vec &x) { return apply_function<double>(std::tan, x); }
  //! Tan function
  inline mat tan(const mat &x) { return apply_function<double>(std::tan, x); }
  //! Inverse sine function
  inline vec asin(const vec &x) { return apply_function<double>(std::asin, x); }
  //! Inverse sine function
  inline mat asin(const mat &x) { return apply_function<double>(std::asin, x); }
  //! Inverse cosine function
  inline vec acos(const vec &x) { return apply_function<double>(std::acos, x); }
  //! Inverse cosine function
  inline mat acos(const mat &x) { return apply_function<double>(std::acos, x); }
  //! Inverse tan function
  inline vec atan(const vec &x) { return apply_function<double>(std::atan, x); }
  //! Inverse tan function
  inline mat atan(const mat &x) { return apply_function<double>(std::atan, x); }
  //! Sinc function, sin(pi*x)/(pi*x)
  inline vec sinc(const vec &x) { return apply_function<double>(sinc, x); }
  //! Sinc function, sin(pi*x)/(pi*x)
  inline mat sinc(const mat &x) { return apply_function<double>(sinc, x); }

  //!@}


  //!\addtogroup hypfunc
  //!@{

  //! Sine hyperbolic function
  inline vec sinh(const vec &x) { return apply_function<double>(std::sinh, x); }
  //! Sine hyperbolic function
  inline mat sinh(const mat &x) { return apply_function<double>(std::sinh, x); }
  //! Cosine hyperbolic function
  inline vec cosh(const vec &x) { return apply_function<double>(std::cosh, x); }
  //! Cosine hyperbolic function
  inline mat cosh(const mat &x) { return apply_function<double>(std::cosh, x); }
  //! Tan hyperbolic function
  inline vec tanh(const vec &x) { return apply_function<double>(std::tanh, x); }
  //! Tan hyperbolic function
  inline mat tanh(const mat &x) { return apply_function<double>(std::tanh, x); }
  //! Inverse sine hyperbolic function
  inline vec asinh(const vec &x) { return apply_function<double>(::asinh, x); }
  //! Inverse sine hyperbolic function
  inline mat asinh(const mat &x) { return apply_function<double>(::asinh, x); }
  //! Inverse cosine hyperbolic function
  inline vec acosh(const vec &x) { return apply_function<double>(::acosh, x); }
  //! Inverse cosine hyperbolic function
  inline mat acosh(const mat &x) { return apply_function<double>(::acosh, x); }
  //! Inverse tan hyperbolic function
  inline vec atanh(const vec &x) { return apply_function<double>(::atanh, x); }
  //! Inverse tan hyperbolic function
  inline mat atanh(const mat &x) { return apply_function<double>(::atanh, x); }

  //!@}

} // namespace itpp

#endif // #ifndef TRIHYPFUNC_H
