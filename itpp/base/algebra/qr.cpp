/*!
 * \file
 * \brief Implementation of QR factorisation functions
 * \author Tony Ottosson, Simon Wood, Adam Piatyszek and Vasek Smidl
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

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#if defined(HAVE_LAPACK)
#  include <itpp/base/algebra/lapack.h>
#endif

#include <itpp/base/algebra/qr.h>
#include <itpp/base/specmat.h>

#include <algorithm>

namespace itpp
{

bool qr(const mat &A, mat &Q, mat &R)
{
#if defined(HAVE_LAPACK)
  int_lapack_t _n, _m, _k, _lwork, info;
  int m = A.rows();
  int n = A.cols();
  int k = std::min(m, n);
  int lwork = n;
  _m = static_cast<int_lapack_t>(m);
  _n = static_cast<int_lapack_t>(n);
  _k = static_cast<int_lapack_t>(k);
  vec tau(k);
  vec work(lwork);

  R = A;

  // perform workspace query for optimum lwork value
  int_lapack_t lwork_tmp = -1;
  dgeqrf_(&_m, &_n, R._data(), &_m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }
  _lwork = static_cast<int_lapack_t>(lwork);
  dgeqrf_(&_m, &_n, R._data(), &_m, tau._data(), work._data(), &_lwork, &info);
  Q = R;
  Q.set_size(m, m, true);

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  // perform workspace query for optimum lwork value
  lwork_tmp = -1;
  dorgqr_(&_m, &_m, &_k, Q._data(), &_m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }
  _lwork = static_cast<int_lapack_t>(lwork);
  dorgqr_(&_m, &_m, &_k, Q._data(), &_m, tau._data(), work._data(), &_lwork,
          &info);

  return (info == 0);
#else
  it_error("LAPACK library is needed to use qr() function");
  return false;
#endif
}

bool qr(const mat &A, mat &R)
{
#if defined(HAVE_LAPACK)
  int_lapack_t _n, _m, _lwork, info;
  int m = A.rows();
  int n = A.cols();
  int k = std::min(m, n);
  int lwork = n;
  _m = static_cast<int_lapack_t>(m);
  _n = static_cast<int_lapack_t>(n);
  vec tau(k);
  vec work(lwork);

  R = A;

  // perform workspace query for optimum lwork value
  int_lapack_t lwork_tmp = -1;
  dgeqrf_(&_m, &_n, R._data(), &_m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }
  _lwork = static_cast<int_lapack_t>(lwork);
  dgeqrf_(&_m, &_n, R._data(), &_m, tau._data(), work._data(), &_lwork, &info);

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  return (info == 0);
#else
  it_error("LAPACK library is needed to use qr() function");
  return false;
#endif
}

bool qr(const mat &A, mat &Q, mat &R, bmat &P)
{
#if defined(HAVE_LAPACK)
  int_lapack_t _n, _m, _k, _lwork, info;
  int m = A.rows();
  int n = A.cols();
  int k = std::min(m, n);
  int lwork = n;
  int_lapack_t *jpvt = new int_lapack_t[n];
  _m = static_cast<int_lapack_t>(m);
  _n = static_cast<int_lapack_t>(n);
  _k = static_cast<int_lapack_t>(k);
  vec tau(k);
  vec work(lwork);

  for (int i = 0; i < n; i++) {
    jpvt[i] = 0;
  }

  R = A;

  // perform workspace query for optimum lwork value
  int_lapack_t lwork_tmp = -1;
  dgeqp3_(&_m, &_n, R._data(), &_m, jpvt, tau._data(), work._data(),
          &lwork_tmp, &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }
  _lwork = static_cast<int_lapack_t>(lwork);
  dgeqp3_(&_m, &_n, R._data(), &_m, jpvt, tau._data(), work._data(),
          &_lwork, &info);
  Q = R;
  Q.set_size(m, m, true);

  // construct permutation matrix
  P = zeros_b(n, n);
  for (int j = 0; j < n; j++)
    P(static_cast<int>(jpvt[j]) - 1, j) = 1;

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  // perform workspace query for optimum lwork value
  lwork_tmp = -1;
  dorgqr_(&_m, &_m, &_k, Q._data(), &_m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }
  _lwork = static_cast<int_lapack_t>(lwork);
  dorgqr_(&_m, &_m, &_k, Q._data(), &_m, tau._data(), work._data(), &_lwork,
          &info);

  delete[] jpvt;

  return (info == 0);
#else
  it_error("LAPACK library is needed to use qr() function");
  return false;
#endif
}

bool qr(const cmat &A, cmat &Q, cmat &R)
{
#if defined(HAVE_LAPACK)
  int_lapack_t _n, _m, _k, _lwork, info;
  int m = A.rows();
  int n = A.cols();
  int k = std::min(m, n);
  int lwork = n;
  _m = static_cast<int_lapack_t>(m);
  _n = static_cast<int_lapack_t>(n);
  _k = static_cast<int_lapack_t>(k);
  cvec tau(k);
  cvec work(lwork);

  R = A;

  // perform workspace query for optimum lwork value
  int_lapack_t lwork_tmp = -1;
  zgeqrf_(&_m, &_n, R._data(), &_m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }
  _lwork = static_cast<int_lapack_t>(lwork);
  zgeqrf_(&_m, &_n, R._data(), &_m, tau._data(), work._data(), &_lwork, &info);

  Q = R;
  Q.set_size(m, m, true);

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  // perform workspace query for optimum lwork value
  lwork_tmp = -1;
  zungqr_(&_m, &_m, &_k, Q._data(), &_m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }
  _lwork = static_cast<int_lapack_t>(lwork);
  zungqr_(&_m, &_m, &_k, Q._data(), &_m, tau._data(), work._data(), &_lwork,
          &info);

  return (info == 0);
#else
  it_error("LAPACK library is needed to use qr() function");
  return false;
#endif
}

bool qr(const cmat &A, cmat &R)
{
#if defined(HAVE_LAPACK)
  int_lapack_t _n, _m, _lwork, info;
  int m = A.rows();
  int n = A.cols();
  int k = std::min(m, n);
  int lwork = n;
  _m = static_cast<int_lapack_t>(m);
  _n = static_cast<int_lapack_t>(n);
  cvec tau(k);
  cvec work(lwork);

  R = A;

  // perform workspace query for optimum lwork value
  int_lapack_t lwork_tmp = -1;
  zgeqrf_(&_m, &_n, R._data(), &_m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }
  _lwork = static_cast<int_lapack_t>(lwork);
  zgeqrf_(&_m, &_n, R._data(), &_m, tau._data(), work._data(), &_lwork, &info);

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  return (info == 0);
#else
  it_error("LAPACK library is needed to use qr() function");
  return false;
#endif
}

bool qr(const cmat &A, cmat &Q, cmat &R, bmat &P)
{
#if defined(HAVE_LAPACK)
  int_lapack_t _n, _m, _k, _lwork, info;
  int m = A.rows();
  int n = A.cols();
  int k = std::min(m, n);
  int lwork = n;
  int_lapack_t *jpvt = new int_lapack_t[n];
  _m = static_cast<int_lapack_t>(m);
  _n = static_cast<int_lapack_t>(n);
  _k = static_cast<int_lapack_t>(k);
  cvec tau(k);
  cvec work(lwork);
  vec rwork(std::max(1, 2*n));

  for (int i = 0; i < n; i++) {
    jpvt[i] = 0;
  }

  R = A;

  // perform workspace query for optimum lwork value
  int_lapack_t lwork_tmp = -1;
  zgeqp3_(&_m, &_n, R._data(), &_m, jpvt, tau._data(), work._data(),
          &lwork_tmp, rwork._data(), &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }
  _lwork = static_cast<int_lapack_t>(lwork);
  zgeqp3_(&_m, &_n, R._data(), &_m, jpvt, tau._data(), work._data(),
          &_lwork, rwork._data(), &info);

  Q = R;
  Q.set_size(m, m, true);

  // construct permutation matrix
  P = zeros_b(n, n);
  for (int j = 0; j < n; j++)
    P(static_cast<int>(jpvt[j]) - 1, j) = 1;

  // construct R
  for (int i = 0; i < m; i++)
    for (int j = 0; j < std::min(i, n); j++)
      R(i, j) = 0;

  // perform workspace query for optimum lwork value
  lwork_tmp = -1;
  zungqr_(&_m, &_m, &_k, Q._data(), &_m, tau._data(), work._data(), &lwork_tmp,
          &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }
  _lwork = static_cast<int_lapack_t>(lwork);
  zungqr_(&_m, &_m, &_k, Q._data(), &_m, tau._data(), work._data(), &_lwork,
          &info);

  delete[] jpvt;

  return (info == 0);
#else
  it_error("LAPACK library is needed to use qr() function");
  return false;
#endif
}

} // namespace itpp
