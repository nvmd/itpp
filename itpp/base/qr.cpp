/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2005 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
 * Tobias Ringström, and Jonas Samuelsson.                                   *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*!
  \file
  \brief Implementation of QR factorisation functions.
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#include "itpp/base/qr.h"
#include "itpp/base/matfunc.h"
#include "itpp/base/lapack.h"

namespace itpp { 

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

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
    it_error("You need to compile IT++ with LAPACK or MKL for qr() to exist");
    return false;
  }

  bool qr(const mat &A, mat &Q, mat &R, bmat &P)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for qr() to exist");
    return false;
  }

  bool qr(const cmat &A, cmat &Q, cmat &R)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for qr() to exist");
    return false;
  }

  bool qr(const cmat &A, cmat &Q, cmat &R, bmat &P)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for qr() to exist");
    return false;
  }


#endif // HAVE_LAPACK or HAVE_MKL



} //namespace itpp
