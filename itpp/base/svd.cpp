/*---------------------------------------------------------------------------*
 *                                   IT++			                         *
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
  \brief Implementation of Singular Value Decompositions
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#include "itpp/base/svd.h"
#include "itpp/base/matfunc.h"
#include "itpp/base/elmatfunc.h"
#include "itpp/base/lapack.h"

namespace itpp { 

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

  bool svd(const mat &A, vec &S)
  {
    char jobu='N', jobvt='N';
    int m, n, lda, ldu, ldvt, lwork, info;
    m = lda = ldu = A.rows();
    n = ldvt = A.cols();
    lwork = std::max(3*std::min(m,n)+std::max(m,n), 5*std::min(m,n));

    mat U, V;
    S.set_size(std::min(m,n), false);
    vec work(lwork);

    mat B(A);

    dgesvd_(&jobu, &jobvt, &m, &n, B._data(), &lda, S._data(), U._data(), &ldu, V._data(), &ldvt, work._data(), &lwork, &info);

    return (info==0);
  }

  bool svd(const cmat &A, vec &S)
  {
    char jobu='N', jobvt='N';
    int m, n, lda, ldu, ldvt, lwork, info;
    m = lda = ldu = A.rows();
    n = ldvt = A.cols();
    lwork = 2*std::min(m,n)+std::max(m,n);

    cvec U, V;
    //U.set_size(m,m, false);
    //V.set_size(n,n, false);
    S.set_size(std::min(m,n), false);
    cvec work(lwork);
    vec rwork(std::max(1, 5*std::min(m, n)));

    cmat B(A);

    zgesvd_(&jobu, &jobvt, &m, &n, B._data(), &lda, S._data(), U._data(), &ldu, V._data(), &ldvt, work._data(), &lwork, rwork._data(), &info);

    return (info==0);
  }

  bool svd(const mat &A, mat &U, vec &S, mat &V)
  {
    char jobu='A', jobvt='A';
    int m, n, lda, ldu, ldvt, lwork, info;
    m = lda = ldu = A.rows();
    n = ldvt = A.cols();
    lwork = std::max(3*std::min(m,n)+std::max(m,n), 5*std::min(m,n));

    U.set_size(m,m, false);
    V.set_size(n,n, false);
    S.set_size(std::min(m,n), false);
    vec work(lwork);

    mat B(A);

    dgesvd_(&jobu, &jobvt, &m, &n, B._data(), &lda, S._data(), U._data(), &ldu, V._data(), &ldvt, work._data(), &lwork, &info);

    V = V.T(); // This is probably slow!!!

    return (info==0);
  }

  bool svd(const cmat &A, cmat &U, vec &S, cmat &V)
  {
    char jobu='A', jobvt='A';
    int m, n, lda, ldu, ldvt, lwork, info;
    m = lda = ldu = A.rows();
    n = ldvt = A.cols();
    lwork = 2*std::min(m,n)+std::max(m,n);

    U.set_size(m,m, false);
    V.set_size(n,n, false);
    S.set_size(std::min(m,n), false);
    cvec work(lwork);
    vec rwork(std::max(1, 5*std::min(m, n)));

    cmat B(A);

    zgesvd_(&jobu, &jobvt, &m, &n, B._data(), &lda, S._data(), U._data(), &ldu, V._data(), &ldvt, work._data(), &lwork, rwork._data(), &info);

    V = V.H(); // This is slow!!!

    return (info==0);
  }

#else
  bool svd(const mat &A, vec &S)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for svd() to exist");
    return false;
  }

  bool svd(const cmat &A, vec &S)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for svd() to exist");
    return false;
  }

  bool svd(const mat &A, mat &U, vec &S, mat &V)
  {   
    it_error("You need to compile IT++ with LAPACK or MKL for svd() to exist");
    return false;
  }

  bool svd(const cmat &A, cmat &U, vec &S, cmat &V)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for svd() to exist");
    return false;   
  }


#endif // HAVE_LAPACK or HAVE_MKL


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
