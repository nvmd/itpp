/*!
 * \file
 * \brief Schur decomposition functions
 * \author Adam Piatyszek
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

#include <itpp/base/algebra/schur.h>


namespace itpp
{

bool schur(const mat &A, mat &U, mat &T)
{
  it_assert_debug(A.rows() == A.cols(), "schur(): Matrix is not square");

#if defined(HAVE_LAPACK)
  int_lapack_t _n, lda, ldvs, sdim, _lwork, info;
  char jobvs = 'V';
  char sort = 'N';
  int n = A.rows();
  int lwork = 3 * n; // This may be choosen better!
  _n = lda = ldvs = static_cast<int_lapack_t>(n);
  _lwork = static_cast<int_lapack_t>(lwork);
  sdim = 0;
  vec wr(n);
  vec wi(n);
  vec work(lwork);

  T.set_size(static_cast<int>(lda), n, false);
  U.set_size(static_cast<int>(ldvs), n, false);

  T = A; // The routine overwrites input matrix with eigenvectors

  dgees_(&jobvs, &sort, 0, &_n, T._data(), &lda, &sdim, wr._data(), wi._data(),
         U._data(), &ldvs, work._data(), &_lwork, 0, &info);

  return (info == 0);
#else
  it_error("LAPACK library is needed to use schur() function");
  return false;
#endif
}


bool schur(const cmat &A, cmat &U, cmat &T)
{
  it_assert_debug(A.rows() == A.cols(), "schur(): Matrix is not square");

#if defined(HAVE_LAPACK)
  int_lapack_t _n, lda, ldvs, sdim, _lwork, info;
  char jobvs = 'V';
  char sort = 'N';
  int n = A.rows();
  int lwork = 2 * n; // This may be choosen better!
  _n = lda = ldvs = static_cast<int_lapack_t>(n);
  _lwork = static_cast<int_lapack_t>(lwork);
  sdim = 0;
  vec rwork(n);
  cvec w(n);
  cvec work(lwork);

  T.set_size(static_cast<int>(lda), n, false);
  U.set_size(static_cast<int>(ldvs), n, false);

  T = A; // The routine overwrites input matrix with eigenvectors

  zgees_(&jobvs, &sort, 0, &_n, T._data(), &lda, &sdim, w._data(), U._data(),
         &ldvs, work._data(), &_lwork, rwork._data(), 0, &info);

  return (info == 0);
#else
  it_error("LAPACK library is needed to use schur() function");
  return false;
#endif
}

mat schur(const mat &A)
{
  mat U, T;
  schur(A, U, T);
  return T;
}


cmat schur(const cmat &A)
{
  cmat U, T;
  schur(A, U, T);
  return T;
}

} // namespace itpp
