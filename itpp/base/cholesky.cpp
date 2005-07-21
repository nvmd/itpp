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
  \brief Implementation of Cholesky factorisation functions
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#include <algorithm>
#include <cassert>
#include <itpp/itconfig.h>
#include <itpp/base/cholesky.h>
#include <itpp/base/lapack.h>

namespace itpp { 

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

  bool chol(const mat &X, mat &F)
  {

    char uplo='U';
    int n, lda, info;
    n = lda = X.rows();

    F = X; // input matrix is overwritten

    dpotrf_(&uplo, &n, F._data(), &lda, &info);

    // Set lower part to zero
    for (int i=0; i<n; i++)
      for(int j=i+1; j<n; j++)
	F(j,i) = 0;

    return (info==0);
  }

  bool chol(const cmat &X, cmat &F)
  {
    char uplo='U';
    int n, lda, info;
    n = lda = X.rows();

    F = X; // input matrix is overwritten

    zpotrf_(&uplo, &n, F._data(), &lda, &info);

    // Set lower part to zero
    for (int i=0; i<n; i++)
      for(int j=i+1; j<n; j++)
	F(j,i) = 0;

    return (info==0);
  }


#else // HAVE_LAPACK or HAVE_MKL
  bool chol(const mat &X, mat &F)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for chol() to exist");
    return false;
  }

  bool chol(const cmat &X, cmat &F)
  {

    it_error("You need to compile IT++ with LAPACK or MKL for chol() to exist");
    return false;
  }

#endif // HAVE_LAPACK or HAVE_MKL


  cmat chol(const cmat &X)
  {
    cmat F;
    if (!chol(X, F)) {
      it_warning("cholesky factorization didn't succeed");
    }

    return F;
  }


  mat chol(const mat &X)
  {
    mat F;
    if (!chol(X, F)) {
      it_warning("cholesky factorization didn't succeed");
    }

    return F;
  }



} //namespace itpp
