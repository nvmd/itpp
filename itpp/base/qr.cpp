/*!
 * \file
 * \brief Implementation of QR factorisation functions
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

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#if defined(HAVE_LAPACK)
#  include <itpp/base/lapack.h>
#endif

#include <itpp/base/qr.h>


namespace itpp { 

#if defined(HAVE_LAPACK)

  bool qr(const mat &A, mat &Q, mat &R)
  {
    int m, n, k, info, lwork, i, j;

    m = A.rows(); n = A.cols();
    lwork = 16*n;
    k = std::min(m,n);
    vec tau(k);
    vec work(lwork);

    R = A;
    dgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork, &info);
    Q = R;

    // construct R
    for (i=0; i<m; i++)
      for (j=0; j<std::min(i,n); j++)
	R(i,j) = 0;

    Q.set_size(m, m, true);
    dorgqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork, &info);

    return (info==0);
  }

  bool qr(const mat &A, mat &Q, mat &R, bmat &P)
  {
    int m, n, k, info, lwork, i, j;

    m = A.rows(); n = A.cols();
    lwork = 16*n;
    k = std::min(m,n);
    vec tau(k);
    vec work(lwork);
    ivec jpvt(n); jpvt.zeros();

    R = A;
    P.set_size(n, n, false); P.zeros();
    dgeqp3_(&m, &n, R._data(), &m, jpvt._data(), tau._data(), work._data(), &lwork, &info);
    Q = R;

    // construct permutation matrix
    for (j=0; j<n; j++)
      P(jpvt(j)-1, j) = 1;

    // construct R
    for (i=0; i<m; i++)
      for (j=0; j<std::min(i,n); j++)
	R(i,j) = 0;

    Q.set_size(m, m, true);
    dorgqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork, &info);

    return (info==0);
  }



  bool qr(const cmat &A, cmat &Q, cmat &R)
  {
    int m, n, k, info, lwork, i, j;

    m = A.rows(); n = A.cols();
    lwork = 16*n;
    k = std::min(m,n);
    cvec tau(k);
    cvec work(lwork);

    R = A;

    zgeqrf_(&m, &n, R._data(), &m, tau._data(), work._data(), &lwork, &info);
    Q = R;

    // construct R
    for (i=0; i<m; i++)
      for (j=0; j<std::min(i,n); j++)
	R(i,j) = 0;


    Q.set_size(m, m, true);
    zungqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork, &info);

    return (info==0);
  }

  bool qr(const cmat &A, cmat &Q, cmat &R, bmat &P)
  {
    int m, n, k, info, lwork, i, j;

    m = A.rows(); n = A.cols();
    lwork = 16*n;
    k = std::min(m,n);
    cvec tau(k);
    cvec work(lwork);
    vec rwork(std::max(1, 2*n));
    ivec jpvt(n); 
    jpvt.zeros();

    R = A;
    P.set_size(n, n, false); P.zeros();

    zgeqp3_(&m, &n, R._data(), &m, jpvt._data(), tau._data(), work._data(), &lwork, rwork._data(), &info);
    Q = R;

    // construct permutation matrix
    for (j=0; j<n; j++)
      P(jpvt(j)-1, j) = 1;

    // construct R
    for (i=0; i<m; i++)
      for (j=0; j<std::min(i,n); j++)
	R(i,j) = 0;


    Q.set_size(m, m, true);
    zungqr_(&m, &m, &k, Q._data(), &m, tau._data(), work._data(), &lwork, &info);

    return (info==0);
  }

#else

  bool qr(const mat &A, mat &Q, mat &R)
  {
    it_error("LAPACK library is needed to use qr() function");
    return false;
  }

  bool qr(const mat &A, mat &Q, mat &R, bmat &P)
  {
    it_error("LAPACK library is needed to use qr() function");
    return false;
  }

  bool qr(const cmat &A, cmat &Q, cmat &R)
  {
    it_error("LAPACK library is needed to use qr() function");
    return false;
  }

  bool qr(const cmat &A, cmat &Q, cmat &R, bmat &P)
  {
    it_error("LAPACK library is needed to use qr() function");
    return false;
  }

#endif // HAVE_LAPACK

} // namespace itpp
