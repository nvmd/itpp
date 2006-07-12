/*!
 * \file
 * \brief Help functions to make functions with vec and mat as arguments
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

#include <itpp/base/help_functions.h>


namespace itpp {

  // Template instantiation of apply_function
  template vec apply_function(double (*f)(double), const vec &v);
  // Template instantiation of apply_function
  template cvec apply_function(std::complex<double> (*f)(const std::complex<double> &),
			       const cvec &v);
  // Template instantiation of apply_function
  template svec apply_function(short (*f)(short), const svec &v);
  // Template instantiation of apply_function
  template ivec apply_function(int (*f)(int), const ivec &v);
  // Template instantiation of apply_function
  template bvec apply_function(bin (*f)(bin), const bvec &v);

  // Template instantiation of apply_function
  template mat apply_function(double (*f)(double), const mat &m);
  // Template instantiation of apply_function
  template cmat apply_function(std::complex<double> (*f)(const std::complex<double> &),
			       const cmat &m);
  // Template instantiation of apply_function
  template smat apply_function(short (*f)(short), const smat &m);
  // Template instantiation of apply_function
  template imat apply_function(int (*f)(int), const imat &m);
  // Template instantiation of apply_function
  template bmat apply_function(bin (*f)(bin), const bmat &m);

  // Template instantiation of apply_function
  template vec apply_function(double (*f)(double, double), const double& x, const vec &v);
  // Template instantiation of apply_function
  template cvec apply_function(std::complex<double> (*f)(const std::complex<double> &, 
							 const std::complex<double> &), 
			       const std::complex<double>& x, const cvec &v);
  // Template instantiation of apply_function
  template svec apply_function(short (*f)(short, short), const short& x, const svec &v);
  // Template instantiation of apply_function
  template ivec apply_function(int (*f)(int, int), const int& x, const ivec &v);
  // Template instantiation of apply_function
  template bvec apply_function(bin (*f)(bin, bin), const bin& x, const bvec &v);

  // Template instantiation of apply_function
  template mat apply_function(double (*f)(double, double), const double& x, const mat &m);
  // Template instantiation of apply_function
  template cmat apply_function(std::complex<double> (*f)(const std::complex<double> &, 
							 const std::complex<double> &), 
			       const std::complex<double>& x, const cmat &m);
  // Template instantiation of apply_function
  template smat apply_function(short (*f)(short, short), const short& x, const smat &m);
  // Template instantiation of apply_function
  template imat apply_function(int (*f)(int, int), const int& x, const imat &m);
  // Template instantiation of apply_function
  template bmat apply_function(bin (*f)(bin, bin), const bin& x, const bmat &m);

  // Template instantiation of apply_function
  template vec apply_function(double (*f)(double, double), const vec &v, const double& x);
  // Template instantiation of apply_function
  template cvec apply_function(std::complex<double> (*f)(const std::complex<double> &,
							 const std::complex<double> &),
			       const cvec &v, const std::complex<double>& x);
  // Template instantiation of apply_function
  template svec apply_function(short (*f)(short, short), const svec &v, const short& x);
  // Template instantiation of apply_function
  template ivec apply_function(int (*f)(int, int), const ivec &v, const int& x);
  // Template instantiation of apply_function
  template bvec apply_function(bin (*f)(bin, bin), const bvec &v, const bin& x);

  // Template instantiation of apply_function
  template mat apply_function(double (*f)(double, double), const mat &m, const double& x);
  // Template instantiation of apply_function
  template cmat apply_function(std::complex<double> (*f)(const std::complex<double> &, 
							 const std::complex<double> &), 
			       const cmat &m, const std::complex<double>& x);
  // Template instantiation of apply_function
  template smat apply_function(short (*f)(short, short), const smat &m, const short& x);
  // Template instantiation of apply_function
  template imat apply_function(int (*f)(int, int), const imat &m, const int& x);
  // Template instantiation of apply_function
  template bmat apply_function(bin (*f)(bin, bin), const bmat &m, const bin& x);

}
