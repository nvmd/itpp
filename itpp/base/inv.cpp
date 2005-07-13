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
  \brief Implementation of matrix inversion routines
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#include "itpp/base/inv.h"
#include "itpp/base/vec.h"
#include "itpp/base/itassert.h"
#include "itpp/base/lapack.h"

namespace itpp { 

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

  bool inv(const mat &X, mat &Y)
  {
    it_assert1(X.rows() == X.cols(), "inv: matrix is not square");

    int m = X.rows(), info, lwork;
    lwork = m; // may be choosen better

    ivec p(m);
    Y = X;
    vec work(lwork);

    dgetrf_(&m, &m, Y._data(), &m, p._data(), &info); // LU-factorization
    if (info!=0)
      return false;

    dgetri_(&m, Y._data(), &m, p._data(), work._data(), &lwork, &info);
    return (info==0);
  }

  bool inv(const cmat &X, cmat &Y)
  {
    it_assert1(X.rows() == X.cols(), "inv: matrix is not square");

    int m = X.rows(), info, lwork;
    lwork = m; // may be choosen better

    ivec p(m);
    Y = X;
    cvec work(lwork);

    zgetrf_(&m, &m, Y._data(), &m, p._data(), &info); // LU-factorization
    if (info!=0)
      return false;

    zgetri_(&m, Y._data(), &m, p._data(), work._data(), &lwork, &info);
    return (info==0);
  }

#else
  bool inv(const mat &X, mat &Y)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for inv() to exist");
    return false;
  }

  bool inv(const cmat &X, cmat &Y)
  {
    it_error("You need to compile IT++ with LAPACK or MKL for inv() to exist");
    return false;
  }

#endif // HAVE_LAPACK or HAVE_MKL


  cmat inv(const cmat &X)
  {
    cmat Y;
    inv(X, Y);
    return Y;
  }


  mat inv(const mat &X)
  {
    mat Y;
    inv(X, Y);
    return Y;
  }

} //namespace itpp
