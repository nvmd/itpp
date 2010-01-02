/*!
 * \file
 * \brief Various functions on vectors and matrices - source file
 * \author Tony Ottosson, Adam Piatyszek, Conrad Sanderson, Mark Dobossy
 *         and Martin Senst
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include <itpp/base/matfunc.h>
#include <itpp/base/algebra/schur.h>
#include <itpp/base/converters.h>
#include <limits>

//! \cond

namespace itpp
{

// Square root of a square matrix (based on Octave sqrtm implementation)
cmat sqrtm(const mat& A)
{
  return sqrtm(to_cmat(A));
}

// Square root of the complex square matrix A
cmat sqrtm(const cmat& A)
{
  cmat U, T;
  schur(A, U, T);

  int n = U.rows();
  cmat R(n, n);

  R.zeros();
  for (int j = 0; j < n; j++)
    R(j, j) = std::sqrt(T(j, j));

  const double fudge = std::sqrt(std::numeric_limits<double>::min());

  for (int p = 0; p < n - 1; p++) {
    for (int i = 0; i < n - (p + 1); i++) {
      const int j = i + p + 1;
      std::complex<double> s = T(i, j);

      for (int k = i + 1; k < j; k++)
        s -= R(i, k) * R(k, j);

      const std::complex<double> d = R(i, i) + R(j, j) + fudge;
      const std::complex<double> conj_d = conj(d);

      R(i, j) = (s * conj_d) / (d * conj_d);
    }
  }

  return U * R * U.H();
}


bool all(const bvec &testvec)
{
  for (int i = 0; i < testvec.length(); i++)
    if (!testvec(i)) return false;
  return true;
}

bool any(const bvec &testvec)
{
  for (int i = 0; i < testvec.length(); i++)
    if (testvec(i)) return true;
  return false;
}

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

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
template cvec cross(const cvec &v1, const cvec &v2);
template ivec cross(const ivec &v1, const ivec &v2);
template svec cross(const svec &v1, const svec &v2);
template bvec cross(const bvec &v1, const bvec &v2);

template vec reverse(const vec &in);
template cvec reverse(const cvec &in);
template svec reverse(const svec &in);
template ivec reverse(const ivec &in);
template bvec reverse(const bvec &in);

template vec zero_pad(const vec &v, int n);
template cvec zero_pad(const cvec &v, int n);
template ivec zero_pad(const ivec &v, int n);
template svec zero_pad(const svec &v, int n);
template bvec zero_pad(const bvec &v, int n);

template vec zero_pad(const vec &v);
template cvec zero_pad(const cvec &v);
template ivec zero_pad(const ivec &v);
template svec zero_pad(const svec &v);
template bvec zero_pad(const bvec &v);

template mat  zero_pad(const mat &, int, int);
template cmat zero_pad(const cmat &, int, int);
template imat zero_pad(const imat &, int, int);
template smat zero_pad(const smat &, int, int);
template bmat zero_pad(const bmat &, int, int);

template vec sum(const mat &m, int dim);
template cvec sum(const cmat &m, int dim);
template svec sum(const smat &m, int dim);
template ivec sum(const imat &m, int dim);
template bvec sum(const bmat &m, int dim);

template double sumsum(const mat &X);
template std::complex<double> sumsum(const cmat &X);
template short sumsum(const smat &X);
template int sumsum(const imat &X);
template bin sumsum(const bmat &X);

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
template cvec prod(const cmat &v, int dim);
template svec prod(const smat &m, int dim);
template ivec prod(const imat &m, int dim);
template bvec prod(const bmat &m, int dim);

template vec diag(const mat &in);
template cvec diag(const cmat &in);

template void diag(const vec &in, mat &m);
template void diag(const cvec &in, cmat &m);

template mat diag(const vec &v, const int K);
template cmat diag(const cvec &v, const int K);

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

template bool is_hermitian(const mat &X);
template bool is_hermitian(const cmat &X);

template bool is_unitary(const mat &X);
template bool is_unitary(const cmat &X);

template vec rvectorize(const mat &m);
template cvec rvectorize(const cmat &m);
template ivec rvectorize(const imat &m);
template svec rvectorize(const smat &m);
template bvec rvectorize(const bmat &m);

template vec cvectorize(const mat &m);
template cvec cvectorize(const cmat &m);
template ivec cvectorize(const imat &m);
template svec cvectorize(const smat &m);
template bvec cvectorize(const bmat &m);

template mat reshape(const mat &m, int rows, int cols);
template cmat reshape(const cmat &m, int rows, int cols);
template imat reshape(const imat &m, int rows, int cols);
template smat reshape(const smat &m, int rows, int cols);
template bmat reshape(const bmat &m, int rows, int cols);

template mat reshape(const vec &m, int rows, int cols);
template cmat reshape(const cvec &m, int rows, int cols);
template imat reshape(const ivec &m, int rows, int cols);
template smat reshape(const svec &m, int rows, int cols);
template bmat reshape(const bvec &m, int rows, int cols);

template mat kron(const mat &X, const mat &Y);
template cmat kron(const cmat &X, const cmat &Y);
template imat kron(const imat &X, const imat &Y);
template smat kron(const smat &X, const smat &Y);
template bmat kron(const bmat &X, const bmat &Y);

template vec repmat(const vec &v, int n);
template cvec repmat(const cvec &v, int n);
template ivec repmat(const ivec &v, int n);
template svec repmat(const svec &v, int n);
template bvec repmat(const bvec &v, int n);

template mat repmat(const vec &v, int m, int n, bool transpose);
template cmat repmat(const cvec &v, int m, int n, bool transpose);
template imat repmat(const ivec &v, int m, int n, bool transpose);
template smat repmat(const svec &v, int m, int n, bool transpose);
template bmat repmat(const bvec &v, int m, int n, bool transpose);

template mat repmat(const mat &data, int m, int n);
template cmat repmat(const cmat &data, int m, int n);
template imat repmat(const imat &data, int m, int n);
template smat repmat(const smat &data, int m, int n);
template bmat repmat(const bmat &data, int m, int n);

} // namespace itpp

//! \endcond
