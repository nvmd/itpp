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
  \brief Implementation of LU factorisation functions.
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#include <itpp/base/lu.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/lapack.h>

namespace itpp { 

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

  bool lu(const mat &X, mat &L, mat &U, ivec &p)
  {
    it_assert1(X.rows() == X.cols(), "lu: matrix is not quadratic");
    //int m, n, lda, info;
    //m = n = lda = X.rows();
    int m = X.rows(), info;

    mat A(X);
    L.set_size(m, m, false);
    U.set_size(m, m, false);
    p.set_size(m, false);

    dgetrf_(&m, &m, A._data(), &m, p._data(), &info);

    for (int i=0; i<m; i++) {
      for (int j=i; j<m; j++) {
	if (i == j) { // diagonal
	  L(i,j) = 1;
	  U(i,j) = A(i,j);
	} else { // upper and lower triangular parts
	  L(i,j) = U(j,i) = 0;
	  L(j,i) = A(j,i);
	  U(i,j) = A(i,j);
	}
      }
    }

    p = p - 1; // Fortran counts from 1

    return (info==0);
  }

  // Slower than not using LAPACK when matrix size smaller than approx 20.
  bool lu(const cmat &X, cmat &L, cmat &U, ivec &p)
  {
    it_assert1(X.rows() == X.cols(), "lu: matrix is not quadratic");
    //int m, n, lda, info;
    //m = n = lda = X.rows();
    int m = X.rows(), info;

    cmat A(X);
    L.set_size(m, m, false);
    U.set_size(m, m, false);
    p.set_size(m, false);

    zgetrf_(&m, &m, A._data(), &m, p._data(), &info);

    for (int i=0; i<m; i++) {
      for (int j=i; j<m; j++) {
	if (i == j) { // diagonal
	  L(i,j) = 1;
	  U(i,j) = A(i,j);
	} else { // upper and lower triangular parts
	  L(i,j) = U(j,i) = 0;
	  L(j,i) = A(j,i);
	  U(i,j) = A(i,j);
	}
      }
    }

    p = p - 1; // Fortran counts from 1

    return (info==0);
  }

#else

  bool lu(const mat &X, mat &L, mat &U, ivec &p)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for lu() to exist");
    return false;   
  }

  bool lu(const cmat &X, cmat &L, cmat &U, ivec &p)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for lu() to exist");
    return false;   
  }

#endif // HAVE_LAPACK or HAVE_MKL


  void interchange_permutations(vec &b, const ivec &p)
  {
    assert(b.size() == p.size());
    double temp;

    for (int k=0; k<b.size(); k++) {
      temp=b(k);
      b(k)=b(p(k));
      b(p(k))=temp;
    }
  }

  bmat permutation_matrix(const ivec &p)
  {
    assert (p.size() > 0);
    int n = p.size(), k;
    bmat P, identity;
    bvec row_k, row_pk;
    identity=eye_b(n);


    for (k=n-1; k>=0; k--) {
      // swap rows k and p(k) in identity
      row_k=identity.get_row(k);
      row_pk=identity.get_row(p(k));
      identity.set_row(k, row_pk);
      identity.set_row(p(k), row_k);

      if (k == n-1) {
	P=identity;
      } else {
	P*=identity;
      }

      // swap back
      identity.set_row(k, row_k);
      identity.set_row(p(k), row_pk);
    }
    return P;
  }

} //namespace itpp
