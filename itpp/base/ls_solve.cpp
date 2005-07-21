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
  \brief Implementation of functions for solving linear equation systems.
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#include <cassert>
#include <itpp/itconfig.h>
#include <itpp/base/ls_solve.h>
#include <itpp/base/lapack.h>



namespace itpp { 


  // ----------- ls_solve_chol -----------------------------------------------------------

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

  bool ls_solve_chol(const mat &A, const vec &b, vec &x)
  {
    int n, lda, ldb, nrhs, info;
    n = lda = ldb = A.rows();
    nrhs = 1;
    char uplo='U';

    it_assert1(A.cols() == n, "ls_solve_chol: System-matrix is not square");
    it_assert1(n == b.size(), "The number of rows in A must equal the length of b!");

    ivec ipiv(n);
    x = b;
    mat Chol = A;

    dposv_(&uplo, &n, &nrhs, Chol._data(), &lda, x._data(), &ldb, &info);

    return (info==0);
  }


  bool ls_solve_chol(const mat &A, const mat &B, mat &X)
  {
    int n, lda, ldb, nrhs, info;
    n = lda = ldb = A.rows();
    nrhs = B.cols();
    char uplo='U';

    it_assert1(A.cols() == n, "ls_solve_chol: System-matrix is not square");
    it_assert1(n == B.rows(), "The number of rows in A must equal the length of B!");

    ivec ipiv(n);
    X = B;
    mat Chol = A;

    dposv_(&uplo, &n, &nrhs, Chol._data(), &lda, X._data(), &ldb, &info);

    return (info==0);
  }

  bool ls_solve_chol(const cmat &A, const cvec &b, cvec &x)
  {
    int n, lda, ldb, nrhs, info;
    n = lda = ldb = A.rows();
    nrhs = 1;
    char uplo='U';

    it_assert1(A.cols() == n, "ls_solve_chol: System-matrix is not square");
    it_assert1(n == b.size(), "The number of rows in A must equal the length of b!");

    ivec ipiv(n);
    x = b;
    cmat Chol = A;

    zposv_(&uplo, &n, &nrhs, Chol._data(), &lda, x._data(), &ldb, &info);

    return (info==0);
  }

  bool ls_solve_chol(const cmat &A, const cmat &B, cmat &X)
  {
    int n, lda, ldb, nrhs, info;
    n = lda = ldb = A.rows();
    nrhs = B.cols();
    char uplo='U';

    it_assert1(A.cols() == n, "ls_solve_chol: System-matrix is not square");
    it_assert1(n == B.rows(), "The number of rows in A must equal the length of B!");

    ivec ipiv(n);
    X = B;
    cmat Chol = A;

    zposv_(&uplo, &n, &nrhs, Chol._data(), &lda, X._data(), &ldb, &info);

    return (info==0);
  }


#else
  bool ls_solve_chol(const mat &A, const vec &b, vec &x)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_chol() to exist");
    return false;
  }

  bool ls_solve_chol(const mat &A, const mat &B, mat &X)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_chol() to exist");
    return false;
  }

  bool ls_solve_chol(const cmat &A, const cvec &b, cvec &x)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_chol() to exist");
    return false;
  }

  bool ls_solve_chol(const cmat &A, const cmat &B, cmat &X)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_chol() to exist");
    return false;
  }

#endif // HAVE_LAPACK or HAVE_MKL


  vec ls_solve_chol(const mat &A, const vec &b)
  {
    vec x;
    bool info = ls_solve_chol(A, b, x);

    it_assert1(info, "ls_solve_chol: Failed solving the system");
    return x;
  }

  mat ls_solve_chol(const mat &A, const mat &B)
  {
    mat X;
    bool info = ls_solve_chol(A, B, X);

    it_assert1(info, "ls_solve_chol: Failed solving the system");
    return X;
  }

  cvec ls_solve_chol(const cmat &A, const cvec &b)
  {
    cvec x;
    bool info = ls_solve_chol(A, b, x);

    it_assert1(info, "ls_solve_chol: Failed solving the system");
    return x;
  }

  cmat ls_solve_chol(const cmat &A, const cmat &B)
  {
    cmat X;
    bool info = ls_solve_chol(A, B, X);

    it_assert1(info, "ls_solve_chol: Failed solving the system");
    return X;
  }




  // --------- ls_solve ---------------------------------------------------------------
#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

  bool ls_solve(const mat &A, const vec &b, vec &x)
  {
    int n, lda, ldb, nrhs, info;
    n = lda = ldb = A.rows();
    nrhs = 1;

    it_assert1(A.cols() == n, "ls_solve: System-matrix is not square");
    it_assert1(n == b.size(), "The number of rows in A must equal the length of b!");

    ivec ipiv(n);
    x = b;
    mat LU = A;

    dgesv_(&n, &nrhs, LU._data(), &lda, ipiv._data(), x._data(), &ldb, &info);

    return (info==0);
  }

  bool ls_solve(const mat &A, const mat &B, mat &X)
  {
    int n, lda, ldb, nrhs, info;
    n = lda = ldb = A.rows();
    nrhs = B.cols();

    it_assert1(A.cols() == n, "ls_solve: System-matrix is not square");
    it_assert1(n == B.rows(), "The number of rows in A must equal the length of B!");

    ivec ipiv(n);
    X = B;
    mat LU = A;

    dgesv_(&n, &nrhs, LU._data(), &lda, ipiv._data(), X._data(), &ldb, &info);

    return (info==0);
  }

  bool ls_solve(const cmat &A, const cvec &b, cvec &x)
  {
    int n, lda, ldb, nrhs, info;
    n = lda = ldb = A.rows();
    nrhs = 1;

    it_assert1(A.cols() == n, "ls_solve: System-matrix is not square");
    it_assert1(n == b.size(), "The number of rows in A must equal the length of b!");

    ivec ipiv(n);
    x = b;
    cmat LU = A;

    zgesv_(&n, &nrhs, LU._data(), &lda, ipiv._data(), x._data(), &ldb, &info);

    return (info==0);
  }

  bool ls_solve(const cmat &A, const cmat &B, cmat &X)
  {
    int n, lda, ldb, nrhs, info;
    n = lda = ldb = A.rows();
    nrhs = B.cols();

    it_assert1(A.cols() == n, "ls_solve: System-matrix is not square");
    it_assert1(n == B.rows(), "The number of rows in A must equal the length of B!");

    ivec ipiv(n);
    X = B;
    cmat LU = A;

    zgesv_(&n, &nrhs, LU._data(), &lda, ipiv._data(), X._data(), &ldb, &info);

    return (info==0);
  }

#else

  bool ls_solve(const mat &A, const vec &b, vec &x)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve() to exist");
    return false;   
  }

  bool ls_solve(const mat &A, const mat &B, mat &X)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve() to exist");
    return false;   
  }

  bool ls_solve(const cmat &A, const cvec &b, cvec &x)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve() to exist");
    return false;   
  }

  bool ls_solve(const cmat &A, const cmat &B, cmat &X)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve() to exist");
    return false;   
  }

#endif // HAVE_LAPACK or HAVE_MKL



  vec ls_solve(const mat &A, const vec &b)
  {
    vec x;
    bool info = ls_solve(A, b, x);

    it_assert1(info, "ls_solve: Failed solving the system");
    return x;
  }

  mat ls_solve(const mat &A, const mat &B)
  {
    mat X;
    bool info = ls_solve(A, B, X);

    it_assert1(info, "ls_solve: Failed solving the system");
    return X;
  }

  cvec ls_solve(const cmat &A, const cvec &b)
  {
    cvec x;
    bool info = ls_solve(A, b, x);

    it_assert1(info, "ls_solve: Failed solving the system");
    return x;
  }

  cmat ls_solve(const cmat &A, const cmat &B)
  {
    cmat X;
    bool info = ls_solve(A, B, X);

    it_assert1(info, "ls_solve: Failed solving the system");
    return X;
  }


  // ----------------- ls_solve_od ------------------------------------------------------------------
#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

  bool ls_solve_od(const mat &A, const vec &b, vec &x)
  {
    int m, n, lda, ldb, nrhs, lwork, info;
    char trans='N';
    m = lda = ldb = A.rows();
    n = A.cols();
    nrhs = 1;
    lwork = n + std::max(m,nrhs);

    it_assert1(m >= n, "The system is under-determined!");
    it_assert1(m == b.size(), "The number of rows in A must equal the length of b!");

    vec work(lwork);
    x = b;
    mat QR = A;

    dgels_(&trans, &m, &n, &nrhs, QR._data(), &lda, x._data(), &ldb, work._data(), &lwork, &info);
    x.set_size(n, true);

    return (info==0);
  }

  bool ls_solve_od(const mat &A, const mat &B, mat &X)
  {
    int m, n, lda, ldb, nrhs, lwork, info;
    char trans='N';
    m = lda = ldb = A.rows();
    n = A.cols();
    nrhs = B.cols();
    lwork = n + std::max(m,nrhs);

    it_assert1(m >= n, "The system is under-determined!");
    it_assert1(m == B.rows(), "The number of rows in A must equal the length of b!");

    vec work(lwork);
    X = B;
    mat QR = A;

    dgels_(&trans, &m, &n, &nrhs, QR._data(), &lda, X._data(), &ldb, work._data(), &lwork, &info);
    X.set_size(n, nrhs, true);

    return (info==0);    
  }

  bool ls_solve_od(const cmat &A, const cvec &b, cvec &x)
  {
    int m, n, lda, ldb, nrhs, lwork, info;
    char trans='N';
    m = lda = ldb = A.rows();
    n = A.cols();
    nrhs = 1;
    lwork = n + std::max(m,nrhs);

    it_assert1(m >= n, "The system is under-determined!");
    it_assert1(m == b.size(), "The number of rows in A must equal the length of b!");

    cvec work(lwork);
    x = b;
    cmat QR = A;

    zgels_(&trans, &m, &n, &nrhs, QR._data(), &lda, x._data(), &ldb, work._data(), &lwork, &info);
    x.set_size(n, true);

    return (info==0);
  }

  bool ls_solve_od(const cmat &A, const cmat &B, cmat &X)
  {
    int m, n, lda, ldb, nrhs, lwork, info;
    char trans='N';
    m = lda = ldb = A.rows();
    n = A.cols();
    nrhs = B.cols();
    lwork = n + std::max(m,nrhs);

    it_assert1(m >= n, "The system is under-determined!");
    it_assert1(m == B.rows(), "The number of rows in A must equal the length of b!");

    cvec work(lwork);
    X = B;
    cmat QR = A;

    zgels_(&trans, &m, &n, &nrhs, QR._data(), &lda, X._data(), &ldb, work._data(), &lwork, &info);
    X.set_size(n, nrhs, true);

    return (info==0);    
  }

#else
  bool ls_solve_od(const mat &A, const vec &b, vec &x)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_od() to exist");
    return false;   
  }

  bool ls_solve_od(const mat &A, const mat &B, mat &X)
  { 
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_od() to exist");
    return false;   
  }

  bool ls_solve_od(const cmat &A, const cvec &b, cvec &x)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_od() to exist");
    return false;   
  }

  bool ls_solve_od(const cmat &A, const cmat &B, cmat &X)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_od() to exist");
    return false;   
  }

#endif // HAVE_LAPACK or HAVE_MKL


  vec ls_solve_od(const mat &A, const vec &b)
  {
    vec x;
    bool info = ls_solve_od(A, b, x);

    it_assert1(info, "ls_solve_od: Failed solving the system");
    return x;
  }

  mat ls_solve_od(const mat &A, const mat &B)
  {
    mat X;
    bool info = ls_solve_od(A, B, X);

    it_assert1(info, "ls_solve_od: Failed solving the system");
    return X;
  }

  cvec ls_solve_od(const cmat &A, const cvec &b)
  {
    cvec x;
    bool info = ls_solve_od(A, b, x);

    it_assert1(info, "ls_solve_od: Failed solving the system");
    return x;
  }

  cmat ls_solve_od(const cmat &A, const cmat &B)
  {
    cmat X;
    bool info = ls_solve_od(A, B, X);

    it_assert1(info, "ls_solve_od: Failed solving the system");
    return X;
  }

  // ------------------- ls_solve_ud -----------------------------------------------------------
#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

  bool ls_solve_ud(const mat &A, const vec &b, vec &x)
  {
    int m, n, lda, ldb, nrhs, lwork, info;
    char trans='N';
    m = lda = A.rows();
    n = A.cols();
    ldb = n;
    nrhs = 1;
    lwork = m + std::max(n,nrhs);

    it_assert1(m < n, "The system is over-determined!");
    it_assert1(m == b.size(), "The number of rows in A must equal the length of b!");

    vec work(lwork);
    x = b;
    x.set_size(n, true);
    mat QR = A;

    dgels_(&trans, &m, &n, &nrhs, QR._data(), &lda, x._data(), &ldb, work._data(), &lwork, &info);

    return (info==0);
  }

  bool ls_solve_ud(const mat &A, const mat &B, mat &X)
  {
    int m, n, lda, ldb, nrhs, lwork, info;
    char trans='N';
    m = lda = A.rows();
    n = A.cols();
    ldb = n;
    nrhs = B.cols();
    lwork = m + std::max(n,nrhs);

    it_assert1(m < n, "The system is over-determined!");
    it_assert1(m == B.rows(), "The number of rows in A must equal the length of b!");

    vec work(lwork);
    X = B;
    X.set_size(n, std::max(m, nrhs), true);
    mat QR = A;

    dgels_(&trans, &m, &n, &nrhs, QR._data(), &lda, X._data(), &ldb, work._data(), &lwork, &info);
    X.set_size(n, nrhs, true);

    return (info==0);    
  }

  bool ls_solve_ud(const cmat &A, const cvec &b, cvec &x)
  {
    int m, n, lda, ldb, nrhs, lwork, info;
    char trans='N';
    m = lda = A.rows();
    n = A.cols();
    ldb = n;
    nrhs = 1;
    lwork = m + std::max(n,nrhs);

    it_assert1(m < n, "The system is over-determined!");
    it_assert1(m == b.size(), "The number of rows in A must equal the length of b!");

    cvec work(lwork);
    x = b;
    x.set_size(n, true);
    cmat QR = A;

    zgels_(&trans, &m, &n, &nrhs, QR._data(), &lda, x._data(), &ldb, work._data(), &lwork, &info);

    return (info==0);
  }

  bool ls_solve_ud(const cmat &A, const cmat &B, cmat &X)
  {
    int m, n, lda, ldb, nrhs, lwork, info;
    char trans='N';
    m = lda = A.rows();
    n = A.cols();
    ldb = n;
    nrhs = B.cols();
    lwork = m + std::max(n,nrhs);

    it_assert1(m < n, "The system is over-determined!");
    it_assert1(m == B.rows(), "The number of rows in A must equal the length of b!");

    cvec work(lwork);
    X = B;
    X.set_size(n, std::max(m, nrhs), true);
    cmat QR = A;

    zgels_(&trans, &m, &n, &nrhs, QR._data(), &lda, X._data(), &ldb, work._data(), &lwork, &info);
    X.set_size(n, nrhs, true);

    return (info==0);    
  }

#else

  bool ls_solve_ud(const mat &A, const vec &b, vec &x)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_ud() to exist");
    return false;   
  }

  bool ls_solve_ud(const mat &A, const mat &B, mat &X)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_ud() to exist");
    return false;   
  }

  bool ls_solve_ud(const cmat &A, const cvec &b, cvec &x)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_ud() to exist");
    return false;   
  }

  bool ls_solve_ud(const cmat &A, const cmat &B, cmat &X)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for ls_solve_ud() to exist");
    return false;   
  }

#endif // HAVE_LAPACK or HAVE_MKL


  vec ls_solve_ud(const mat &A, const vec &b)
  {
    vec x;
    bool info = ls_solve_ud(A, b, x);

    it_assert1(info, "ls_solve_ud: Failed solving the system");
    return x;
  }

  mat ls_solve_ud(const mat &A, const mat &B)
  {
    mat X;
    bool info = ls_solve_ud(A, B, X);

    it_assert1(info, "ls_solve_ud: Failed solving the system");
    return X;
  }

  cvec ls_solve_ud(const cmat &A, const cvec &b)
  {
    cvec x;
    bool info = ls_solve_ud(A, b, x);

    it_assert1(info, "ls_solve_ud: Failed solving the system");
    return x;
  }

  cmat ls_solve_ud(const cmat &A, const cmat &B)
  {
    cmat X;
    bool info = ls_solve_ud(A, B, X);

    it_assert1(info, "ls_solve_ud: Failed solving the system");
    return X;
  }


  // ---------------------- backslash ----------------------------------------------------------
  bool backslash(const mat &A, const vec &b, vec &x)
  {
    int m=A.rows(), n=A.cols();
    bool info;

    if (m == n)
      info = ls_solve(A,b,x);
    else if (m > n)
      info = ls_solve_od(A,b,x);
    else
      info = ls_solve_ud(A,b,x);
    
    return info;
  }

  vec backslash(const mat &A, const vec &b)
  {
    vec x;
    it_assert1( backslash(A, b, x), "backslash: solution was not found");

    return x;
  }


  bool backslash(const mat &A, const mat &B, mat &X)
  {
    int m=A.rows(), n=A.cols();
    bool info;

    if (m == n)
      info = ls_solve(A, B, X);
    else if (m > n)
      info = ls_solve_od(A, B, X);
    else
      info = ls_solve_ud(A, B, X);
    
    return info;
  }

  mat backslash(const mat &A, const mat &B)
  {
    mat X;
    it_assert1( backslash(A, B, X), "backslash: solution was not found");

    return X;
  }



  bool backslash(const cmat &A, const cvec &b, cvec &x)
  {
    int m=A.rows(), n=A.cols();
    bool info;

    if (m == n)
      info = ls_solve(A,b,x);
    else if (m > n)
      info = ls_solve_od(A,b,x);
    else
      info = ls_solve_ud(A,b,x);
    
    return info;
  }

  cvec backslash(const cmat &A, const cvec &b)
  {
    cvec x;
    it_assert1( backslash(A, b, x), "backslash: solution was not found");

    return x;
  }


  bool backslash(const cmat &A, const cmat &B, cmat &X)
  {
    int m=A.rows(), n=A.cols();
    bool info;

    if (m == n)
      info = ls_solve(A, B, X);
    else if (m > n)
      info = ls_solve_od(A, B, X);
    else
      info = ls_solve_ud(A, B, X);
    
    return info;
  }

  cmat backslash(const cmat &A, const cmat &B)
  {
    cmat X;
    it_assert1( backslash(A, B, X), "backslash: solution was not found");

    return X;
  }


  // ----------------------------------------------------------------------------------------------
  vec forward_substitution(const mat &L, const vec &b)
  {
    int n = L.rows();
    vec x(n);
  
    forward_substitution(L, b, x);

    return x;
  }

  void forward_substitution(const mat &L, const vec &b, vec &x)
  {
    assert( L.rows() == L.cols() && L.cols() == b.size() && b.size() == x.size() );
    int n = L.rows(), i, j;
    double temp;

    x(0)=b(0)/L(0,0);
    for (i=1;i<n;i++) {
      // Should be: x(i)=((b(i)-L(i,i,0,i-1)*x(0,i-1))/L(i,i))(0); but this is to slow.
      //i_pos=i*L._row_offset();
      temp=0;
      for (j=0; j<i; j++) {
	temp += L._elem(i,j) * x(j);
	//temp+=L._data()[i_pos+j]*x(j);
      }
      x(i) = (b(i)-temp)/L._elem(i,i);
      //x(i)=(b(i)-temp)/L._data()[i_pos+i];
    }
  }

  vec forward_substitution(const mat &L, int p, const vec &b)
  {
    int n = L.rows();
    vec x(n);
  
    forward_substitution(L, p, b, x);

    return x;
  }

  void forward_substitution(const mat &L, int p, const vec &b, vec &x)
  {
    assert( L.rows() == L.cols() && L.cols() == b.size() && b.size() == x.size() && p <= L.rows()/2 );
    int n = L.rows(), i, j;

    x=b;
  
    for (j=0;j<n;j++) {
      x(j)/=L(j,j);
      for (i=j+1;i<std::min(j+p+1,n);i++) {
	x(i)-=L(i,j)*x(j);
      }
    }
  }

  vec backward_substitution(const mat &U, const vec &b)
  {
    vec x(U.rows());
    backward_substitution(U, b, x);

    return x;
  }

  void backward_substitution(const mat &U, const vec &b, vec &x)
  {
    assert( U.rows() == U.cols() && U.cols() == b.size() && b.size() == x.size() );
    int n = U.rows(), i, j;
    double temp;

    x(n-1)=b(n-1)/U(n-1,n-1);
    for (i=n-2; i>=0; i--) {
      // Should be: x(i)=((b(i)-U(i,i,i+1,n-1)*x(i+1,n-1))/U(i,i))(0); but this is too slow.
      temp=0;
      //i_pos=i*U._row_offset();
      for (j=i+1; j<n; j++) {
	temp += U._elem(i,j) * x(j);
	//temp+=U._data()[i_pos+j]*x(j);
      }
      x(i) = (b(i)-temp)/U._elem(i,i);
      //x(i)=(b(i)-temp)/U._data()[i_pos+i];
    }
  }

  vec backward_substitution(const mat &U, int q, const vec &b)
  {
    vec x(U.rows());
    backward_substitution(U, q, b, x);

    return x;
  }

  void backward_substitution(const mat &U, int q, const vec &b, vec &x)
  {
    assert( U.rows() == U.cols() && U.cols() == b.size() && b.size() == x.size() && q <= U.rows()/2);
    int n = U.rows(), i, j;

    x=b;
  
    for (j=n-1; j>=0; j--) {
      x(j) /= U(j,j);
      for (i=std::max(0,j-q); i<j; i++) {
	x(i)-=U(i,j)*x(j);
      }
    }
  }

} //namespace itpp
