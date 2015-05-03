/*!
 * \file
 * \brief Implementation of Singular Value Decompositions
 * \author Tony Ottosson, Simon Wood and Adam Piatyszek
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

#include <itpp/base/algebra/svd.h>
#include <algorithm>

namespace itpp
{

bool svd(const mat &A, vec &S)
{
#if defined(HAVE_LAPACK)
  int_lapack_t _n, _m, lda, ldu, ldvt, _lwork, info;
  char jobu = 'N';
  char jobvt = 'N';
  int m = A.rows();
  int n = A.cols();
  int lwork = std::max(3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n));
  _m = lda = ldu = static_cast<int_lapack_t>(m);
  _n = ldvt = static_cast<int_lapack_t>(n);

  mat U;
  mat V;
  S.set_size(std::min(m, n), false);
  vec work(lwork);

  mat B(A);

  // The theoretical calculation of lwork above results in the minimum size
  // needed for dgesvd_ to run to completion without having memory errors.
  // For speed improvement it is best to set lwork=-1 and have dgesvd_
  // calculate the best workspace requirement.
  int_lapack_t lwork_tmp = -1;
  dgesvd_(&jobu, &jobvt, &_m, &_n, B._data(), &lda, S._data(), U._data(), &ldu,
          V._data(), &ldvt, work._data(), &lwork_tmp, &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }

  _lwork = static_cast<int_lapack_t>(lwork);
  dgesvd_(&jobu, &jobvt, &_m, &_n, B._data(), &lda, S._data(), U._data(), &ldu,
          V._data(), &ldvt, work._data(), &_lwork, &info);

  return (info == 0);
#else
  it_error("LAPACK library is needed to use svd() function");
  return false;
#endif
}

bool svd(const cmat &A, vec &S)
{
#if defined(HAVE_LAPACK)
  int_lapack_t _n, _m, lda, ldu, ldvt, _lwork, info;
  char jobu = 'N';
  char jobvt = 'N';
  int m = A.rows();
  int n = A.cols();
  int lwork = 2 * std::min(m, n) + std::max(m, n);
  _m = lda = ldu = static_cast<int_lapack_t>(m);
  _n = ldvt = static_cast<int_lapack_t>(n);

  cvec U;
  cvec V;
  S.set_size(std::min(m, n), false);
  cvec work(lwork);
  vec rwork(5*std::min(m, n));

  cmat B(A);

  // The theoretical calculation of lwork above results in the minimum size
  // needed for zgesvd_ to run to completion without having memory errors.
  // For speed improvement it is best to set lwork=-1 and have zgesvd_
  // calculate the best workspace requirement.
  int_lapack_t lwork_tmp = -1;
  zgesvd_(&jobu, &jobvt, &_m, &_n, B._data(), &lda, S._data(), U._data(), &ldu,
          V._data(), &ldvt, work._data(), &lwork_tmp, rwork._data(), &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }

  _lwork = static_cast<int_lapack_t>(lwork);
  zgesvd_(&jobu, &jobvt, &_m, &_n, B._data(), &lda, S._data(), U._data(), &ldu,
          V._data(), &ldvt, work._data(), &_lwork, rwork._data(), &info);

  return (info == 0);
#else
  it_error("LAPACK library is needed to use svd() function");
  return false;
#endif
}

bool svd(const mat &A, mat &U, vec &S, mat &V)
{
#if defined(HAVE_LAPACK)
  int_lapack_t _n, _m, lda, ldu, ldvt, _lwork, info;
  char jobu = 'A';
  char jobvt = 'A';
  int m = A.rows();
  int n = A.cols();
  int lwork = std::max(3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n));
  _m = lda = ldu = static_cast<int_lapack_t>(m);
  _n = ldvt = static_cast<int_lapack_t>(n);

  U.set_size(m, m, false);
  V.set_size(n, n, false);
  S.set_size(std::min(m, n), false);
  vec work(lwork);

  mat B(A);

  // The theoretical calculation of lwork above results in the minimum size
  // needed for dgesvd_ to run to completion without having memory errors.
  // For speed improvement it is best to set lwork=-1 and have dgesvd_
  // calculate the best workspace requirement.
  int_lapack_t lwork_tmp = -1;
  dgesvd_(&jobu, &jobvt, &_m, &_n, B._data(), &lda, S._data(), U._data(), &ldu,
          V._data(), &ldvt, work._data(), &lwork_tmp, &info);
  if (info == 0) {
    lwork = static_cast<int>(work(0));
    work.set_size(lwork, false);
  }

  _lwork = static_cast<int_lapack_t>(lwork);
  dgesvd_(&jobu, &jobvt, &_m, &_n, B._data(), &lda, S._data(), U._data(), &ldu,
          V._data(), &ldvt, work._data(), &_lwork, &info);

  V = V.T(); // This is probably slow!!!

  return (info == 0);
#else
  it_error("LAPACK library is needed to use svd() function");
  return false;
#endif
}

bool svd(const cmat &A, cmat &U, vec &S, cmat &V)
{
#if defined(HAVE_LAPACK)
  int_lapack_t _n, _m, lda, ldu, ldvt, _lwork, info;
  char jobu = 'A';
  char jobvt = 'A';
  int m = A.rows();
  int n = A.cols();
  int lwork = 2 * std::min(m, n) + std::max(m, n);
  _m = lda = ldu = static_cast<int_lapack_t>(m);
  _n = ldvt = static_cast<int_lapack_t>(n);

  U.set_size(m, m, false);
  V.set_size(n, n, false);
  S.set_size(std::min(m, n), false);
  cvec work(lwork);
  vec rwork(5 * std::min(m, n));

  cmat B(A);

  // The theoretical calculation of lwork above results in the minimum size
  // needed for zgesvd_ to run to completion without having memory errors.
  // For speed improvement it is best to set lwork=-1 and have zgesvd_
  // calculate the best workspace requirement.
  int_lapack_t lwork_tmp = -1;
  zgesvd_(&jobu, &jobvt, &_m, &_n, B._data(), &lda, S._data(), U._data(), &ldu,
          V._data(), &ldvt, work._data(), &lwork_tmp, rwork._data(), &info);
  if (info == 0) {
    lwork = static_cast<int>(real(work(0)));
    work.set_size(lwork, false);
  }

  _lwork = static_cast<int_lapack_t>(lwork);
  zgesvd_(&jobu, &jobvt, &_m, &_n, B._data(), &lda, S._data(), U._data(), &ldu,
          V._data(), &ldvt, work._data(), &_lwork, rwork._data(), &info);

  V = V.H(); // This is slow!!!

  return (info == 0);
#else
  it_error("LAPACK library is needed to use svd() function");
  return false;
#endif
}

vec svd(const mat &A)
{
  vec S;
  svd(A, S);
  return S;
}

vec svd(const cmat &A)
{
  vec S;
  svd(A, S);
  return S;
}

} //namespace itpp
