/*!
 * \file
 * \brief Implementation of functions on vectors and matrices
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

#include <itpp/base/matfunc.h>

namespace itpp {

  // ---------------------- Instantiations ------------------------------------

  template int length(const vec &v);
  template int length(const cvec &v);
  template int length(const svec &v);
  template int length(const ivec &v);
  template int length(const bvec &v);


  template double sum(const vec &v);
  template std::complex<double> sum(const cvec &v);
  template short sum(const svec &v);
  template int sum(const ivec &v);
  template bin sum(const bvec &v);

  template double sum_sqr(const vec &v);
  template std::complex<double> sum_sqr(const cvec &v);
  template short sum_sqr(const svec &v);
  template int sum_sqr(const ivec &v);
  template bin sum_sqr(const bvec &v);

  template vec cumsum(const vec &v);
  template cvec cumsum(const cvec &v);
  template svec cumsum(const svec &v);
  template ivec cumsum(const ivec &v);
  template bvec cumsum(const bvec &v);

  template double prod(const vec &v);
  template std::complex<double> prod(const cvec &v);
  template short prod(const svec &v);
  template int prod(const ivec &v);
  template bin prod(const bvec &v);

  template vec cross(const vec &v1, const vec &v2);
  template ivec cross(const ivec &v1, const ivec &v2);
  template svec cross(const svec &v1, const svec &v2);

  template vec reverse(const vec &in);
  template cvec reverse(const cvec &in);
  template svec reverse(const svec &in);
  template ivec reverse(const ivec &in);
  template bvec reverse(const bvec &in);

  template vec repeat(const vec &v, int norepeats);
  template cvec repeat(const cvec &v, int norepeats);
  template svec repeat(const svec &v, int norepeats);
  template ivec repeat(const ivec &v, int norepeats);
  template bvec repeat(const bvec &v, int norepeats);

  template vec apply_function(float (*f)(float), const vec &data);
  template vec apply_function(double (*f)(double), const vec &data);
  template cvec apply_function(std::complex<double> (*f)(std::complex<double>), const cvec &data);
  template svec apply_function(short (*f)(short), const svec &data);
  template ivec apply_function(int (*f)(int), const ivec &data);
  template bvec apply_function(bin (*f)(bin), const bvec &data);


  template ivec zero_pad(const ivec &v, int n);
  template vec zero_pad(const vec &v, int n);
  template cvec zero_pad(const cvec &v, int n);
  template bvec zero_pad(const bvec &v, int n);

  template ivec zero_pad(const ivec &v);
  template vec zero_pad(const vec &v);
  template cvec zero_pad(const cvec &v);
  template bvec zero_pad(const bvec &v);

  template mat  zero_pad(const mat &, int, int);
  template cmat zero_pad(const cmat &, int, int);
  template imat zero_pad(const imat &, int, int);
  template bmat zero_pad(const bmat &, int, int);

  template vec sum(const mat &m, int dim);
  template cvec sum(const cmat &m, int dim);
  template svec sum(const smat &m, int dim);
  template ivec sum(const imat &m, int dim);
  template bvec sum(const bmat &m, int dim);

  template vec sum_sqr(const mat & m, int dim);
  template cvec sum_sqr(const cmat &m, int dim);
  template svec sum_sqr(const smat &m, int dim);
  template ivec sum_sqr(const imat &m, int dim);
  template bvec sum_sqr(const bmat &m, int dim);

  template mat cumsum(const mat &m, int dim);
  template cmat cumsum(const cmat &m, int dim);
  template smat cumsum(const smat &m, int dim);
  template imat cumsum(const imat &m, int dim);
  template bmat cumsum(const bmat &m, int dim);

  template vec prod(const mat &m, int dim);
  // Template instantiation of product
  template cvec prod(const cmat &v, int dim);
  template svec prod(const smat &m, int dim);
  template ivec prod(const imat &m, int dim);

  template vec diag(const mat &in);
  template cvec diag(const cmat &in);

  template void diag(const vec &in, mat &m);
  template void diag(const cvec &in, cmat &m);

  template mat diag(const vec &v);
  template cmat diag(const cvec &v);

  template mat bidiag(const vec &, const vec &);
  template cmat bidiag(const cvec &, const cvec &);

  template void bidiag(const vec &, const vec &, mat &);
  template void bidiag(const cvec &, const cvec &, cmat &);

  template void bidiag(const mat &, vec &, vec &);
  template void bidiag(const cmat &, cvec &, cvec &);

  template mat tridiag(const vec &main, const vec &, const vec &);
  template cmat tridiag(const cvec &main, const cvec &, const cvec &);

  template void tridiag(const vec &main, const vec &, const vec &, mat &);
  template void tridiag(const cvec &main, const cvec &, const cvec &, cmat &);

  template void tridiag(const mat &m, vec &, vec &, vec &);
  template void tridiag(const cmat &m, cvec &, cvec &, cvec &);

  template double trace(const mat &in);
  template std::complex<double> trace(const cmat &in);
  template short trace(const smat &in);
  template int trace(const imat &in);
  template bin trace(const bmat &in);

  template void transpose(const mat &m, mat &out);
  template void transpose(const cmat &m, cmat &out);
  template void transpose(const smat &m, smat &out);
  template void transpose(const imat &m, imat &out);
  template void transpose(const bmat &m, bmat &out);

  template mat transpose(const mat &m);
  template cmat transpose(const cmat &m);
  template smat transpose(const smat &m);
  template imat transpose(const imat &m);
  template bmat transpose(const bmat &m);


  template void hermitian_transpose(const mat &m, mat &out);
  template void hermitian_transpose(const cmat &m, cmat &out);
  template void hermitian_transpose(const smat &m, smat &out);
  template void hermitian_transpose(const imat &m, imat &out);
  template void hermitian_transpose(const bmat &m, bmat &out);

  template mat hermitian_transpose(const mat &m);
  template cmat hermitian_transpose(const cmat &m);
  template smat hermitian_transpose(const smat &m);
  template imat hermitian_transpose(const imat &m);
  template bmat hermitian_transpose(const bmat &m);

  template mat repeat(const mat &m, int norepeats);
  template cmat repeat(const cmat &m, int norepeats);
  template smat repeat(const smat &m, int norepeats);
  template imat repeat(const imat &m, int norepeats);
  template bmat repeat(const bmat &m, int norepeats);

  template mat apply_function(float (*f)(float), const mat &data);
  template mat apply_function(double (*f)(double), const mat &data);
  template cmat apply_function(std::complex<double> (*f)(std::complex<double>), const cmat &data);
  template smat apply_function(short (*f)(short), const smat &data);
  template imat apply_function(int (*f)(int), const imat &data);
  template bmat apply_function(bin (*f)(bin), const bmat &data);

  template  vec rvectorize(const  mat &m);
  template cvec rvectorize(const cmat &m);
  template  ivec rvectorize(const  imat &m);
  template  bvec rvectorize(const  bmat &m);

  template  vec cvectorize(const  mat &m);
  template cvec cvectorize(const cmat &m);
  template  ivec cvectorize(const  imat &m);
  template  bvec cvectorize(const  bmat &m);

  template  mat reshape(const  mat &m, int rows, int cols);
  template cmat reshape(const cmat &m, int rows, int cols);
  template  imat reshape(const  imat &m, int rows, int cols);
  template  bmat reshape(const  bmat &m, int rows, int cols);

  template  mat reshape(const  vec &m, int rows, int cols);
  template cmat reshape(const cvec &m, int rows, int cols);
  template  imat reshape(const  ivec &m, int rows, int cols);
  template  bmat reshape(const  bvec &m, int rows, int cols);

  template vec upsample(const vec &v, int usf);
  template cvec upsample(const cvec &v, int usf);
  template svec upsample(const svec &v, int usf);
  template ivec upsample(const ivec &v, int usf);
  template bvec upsample(const bvec &v, int usf);

  template mat upsample(const mat &v, int usf);
  template cmat upsample(const cmat &v, int usf);
  template smat upsample(const smat &v, int usf);
  template imat upsample(const imat &v, int usf);
  template bmat upsample(const bmat &v, int usf);

  template void upsample(const vec &v, int usf,  vec & u);
  template void upsample(const cvec &v, int usf,  cvec & u);
  template void upsample(const svec &v, int usf,  svec & u);
  template void upsample(const ivec &v, int usf,  ivec & u);
  template void upsample(const bvec &v, int usf,  bvec & u);

  template void upsample(const mat &v, int usf,  mat & u);
  template void upsample(const cmat &v, int usf,  cmat & u);
  template void upsample(const smat &v, int usf,  smat & u);
  template void upsample(const imat &v, int usf,  imat & u);
  template void upsample(const bmat &v, int usf,  bmat & u);
 
  template vec lininterp(const vec &v, int usf);
  template cvec lininterp(const cvec &v, int usf);

  template mat lininterp(const mat &v, int usf);
  template cmat lininterp(const cmat &v, int usf);

  template void lininterp(const vec &v, int usf,  vec & u);
  template void lininterp(const cvec &v, int usf,  cvec & u);

  template void lininterp(const mat &v, int usf,  mat & u);
  template void lininterp(const cmat &v, int usf,  cmat & u);

  template mat lininterp(const mat &m, const double f_base, const double f_ups, const int nrof_samples, const double t_start);
  template cmat lininterp(const cmat &m, const double f_base, const double f_ups, const int nrof_samples, const double t_start);

  template vec lininterp(const vec &v, const double f_base, const double f_ups, const int nrof_samples, const double t_start);
  template cvec lininterp(const cvec &v, const double f_base, const double f_ups, const int nrof_samples, const double t_start);

} // namespace itpp
