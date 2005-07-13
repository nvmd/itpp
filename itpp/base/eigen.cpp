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
  \brief Eigenvalue decomposition functions.
  \author Tony Ottosson

  1.11

  2003/06/19 09:03:05
*/

#include "itpp/base/eigen.h"
#include "itpp/base/matfunc.h"
#include "itpp/base/lapack.h"

namespace itpp { 

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

  bool eig_sym(const mat &A, vec &d, mat &V)
  {
    it_assert1(A.rows() == A.cols(), "eig_sym: Matrix is not symmetric");

    // Test for symmetric?

    char jobz='V', uplo='U';
    int n, lda, lwork, info;
    n = lda = A.rows();
    lwork = std::max(1,3*n-1); // This may be choosen better!

    d.set_size(n, false);
    vec work(lwork);

    V = A; // The routine overwrites input matrix with eigenvectors

    dsyev_(&jobz, &uplo, &n, V._data(), &lda, d._data(), work._data(), &lwork, &info);

    return (info==0);
  }

  bool eig_sym(const mat &A, vec &d)
  {
    it_assert1(A.rows() == A.cols(), "eig_sym: Matrix is not symmetric");

    // Test for symmetric?

    char jobz='N', uplo='U';
    int n, lda, lwork, info;
    n = lda = A.rows();
    lwork = std::max(1,3*n-1); // This may be choosen better!

    d.set_size(n, false);
    vec work(lwork);

    mat B(A); // The routine overwrites input matrix

    dsyev_(&jobz, &uplo, &n, B._data(), &lda, d._data(), work._data(), &lwork, &info);

    return (info==0);
  }

  bool eig_sym(const cmat &A, vec &d, cmat &V)
  {
    it_assert1(A.rows() == A.cols(), "eig_sym: Matrix is not hermitian");

    // Test for symmetric?

    char jobz='V', uplo='U';
    int n, lda, lwork, info;
    n = lda = A.rows();
    lwork = std::max(1,2*n-1); // This may be choosen better!

    d.set_size(n, false);
    cvec work(lwork);
    vec rwork(std::max(1,3*n-2)); // This may be choosen better!

    V = A; // The routine overwrites input matrix with eigenvectors

    zheev_(&jobz, &uplo, &n, V._data(), &lda, d._data(), work._data(), &lwork, rwork._data(), &info);

    return (info==0);
  }

  bool eig_sym(const cmat &A, vec &d)
  {
    it_assert1(A.rows() == A.cols(), "eig_sym: Matrix is not hermitian");

    // Test for symmetric?

    char jobz='N', uplo='U';
    int n, lda, lwork, info;
    n = lda = A.rows();
    lwork = std::max(1,2*n-1); // This may be choosen better!

    d.set_size(n, false);
    cvec work(lwork);
    vec rwork(std::max(1,3*n-2)); // This may be choosen better!

    cmat B(A); // The routine overwrites input matrix

    zheev_(&jobz, &uplo, &n, B._data(), &lda, d._data(), work._data(), &lwork, rwork._data(), &info);

    return (info==0);
  }


  // Non-symmetric matrix
  bool eig(const mat &A, cvec &d, cmat &V)
  {
    it_assert1(A.rows() == A.cols(), "eig: Matrix is not square");

    char jobvl='N', jobvr='V';
    int n, lda, ldvl, ldvr, lwork, info;
    n = lda = A.rows();
    ldvl = 1; ldvr = n;
    lwork = std::max(1,4*n); // This may be choosen better!

    vec work(lwork);
    vec rwork(std::max(1,2*n)); // This may be choosen better
    vec wr(n), wi(n);
    mat vl, vr(n,n);

    mat B(A); // The routine overwrites input matrix

    dgeev_(&jobvl, &jobvr, &n, B._data(), &lda, wr._data(), wi._data(), vl._data(), &ldvl, vr._data(), &ldvr, work._data(), &lwork, &info);

    d = to_cvec(wr, wi);

    // Fix V
    V.set_size(n, n, false);
    for (int j=0; j<n; j++) {
      // if d(j) and d(j+1) are complex conjugate pairs, treat special
		if( (j<n-1) && d(j) == std::conj(d(j+1))) {
	V.set_col(j, to_cvec(vr.get_col(j), vr.get_col(j+1)) );
	V.set_col(j+1, to_cvec(vr.get_col(j), -vr.get_col(j+1)) );
	j++;
      } else {
	V.set_col(j, to_cvec(vr.get_col(j)) );
      }
    }

    return (info==0);
  }

  // Non-symmetric matrix
  bool eig(const mat &A, cvec &d)
  {
    it_assert1(A.rows() == A.cols(), "eig: Matrix is not square");

    char jobvl='N', jobvr='N';
    int n, lda, ldvl, ldvr, lwork, info;
    n = lda = A.rows();
    ldvl = 1; ldvr = 1;
    lwork = std::max(1,4*n); // This may be choosen better!

    vec work(lwork);
    vec rwork(std::max(1,2*n)); // This may be choosen better
    vec wr(n), wi(n);
    mat vl, vr;

    mat B(A); // The routine overwrites input matrix

    dgeev_(&jobvl, &jobvr, &n, B._data(), &lda, wr._data(), wi._data(), vl._data(), &ldvl, vr._data(), &ldvr, work._data(), &lwork, &info);

    d = to_cvec(wr, wi);

    return (info==0);
  }

  bool eig(const cmat &A, cvec &d, cmat &V)
  {
    it_assert1(A.rows() == A.cols(), "eig: Matrix is not square");

    char jobvl='N', jobvr='V';
    int n, lda, ldvl, ldvr, lwork, info;
    n = lda = A.rows();
    ldvl = 1; ldvr = n;
    lwork = std::max(1,2*n); // This may be choosen better!

    d.set_size(n, false);
    V.set_size(n, n, false);
    cvec work(lwork);
    vec rwork(std::max(1,2*n)); // This may be choosen better!
    cmat vl;

    cmat B(A); // The routine overwrites input matrix

    zgeev_(&jobvl, &jobvr, &n, B._data(), &lda, d._data(), vl._data(), &ldvl, V._data(), &ldvr, work._data(), &lwork, rwork._data(), &info);


    return (info==0);
  }

  bool eig(const cmat &A, cvec &d)
  {
    it_assert1(A.rows() == A.cols(), "eig: Matrix is not square");

    char jobvl='N', jobvr='N';
    int n, lda, ldvl, ldvr, lwork, info;
    n = lda = A.rows();
    ldvl = 1; ldvr = 1;
    lwork = std::max(1,2*n); // This may be choosen better!

    d.set_size(n, false);
    cvec work(lwork);
    vec rwork(std::max(1,2*n)); // This may be choosen better!
    cmat vl, vr;

    cmat B(A); // The routine overwrites input matrix

    zgeev_(&jobvl, &jobvr, &n, B._data(), &lda, d._data(), vl._data(), &ldvl, vr._data(), &ldvr, work._data(), &lwork, rwork._data(), &info);


    return (info==0);
  }


#else
  bool eig_sym(const mat &A, vec &d, mat &V)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for eig_sym() to exist");
    return false;   
  }

  bool eig_sym(const mat &A, vec &d)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for eig_sym() to exist");
    return false;   
  }

  bool eig_sym(const cmat &A, vec &d, cmat &V)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for eig_sym() to exist");
    return false;   
  }

  bool eig_sym(const cmat &A, vec &d)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for eig_sym() to exist");
    return false;   
  }


  bool eig(const mat &A, cvec &d, cmat &V)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for eig() to exist");
    return false;   
  }

  bool eig(const mat &A, cvec &d)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for eig() to exist");
    return false;   
  }

  bool eig(const cmat &A, cvec &d, cmat &V)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for eig() to exist");
    return false;   
  }

  bool eig(const cmat &A, cvec &d)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for eig() to exist");
    return false;   
  }


#endif // HAVE_LAPACK or HAVE_MKL




  vec eig_sym(const mat &A)
  {
    vec d;
    eig_sym(A, d);
    return d;
  }

  vec eig_sym(const cmat &A)
  {
    vec d;
    eig_sym(A, d);
    return d;
  }

  cvec eig(const mat &A)
  {
    cvec d;
    eig(A, d);
    return d;
  }

  cvec eig(const cmat &A)
  {
    cvec d;
    eig(A, d);
    return d;
  }


} //namespace itpp
