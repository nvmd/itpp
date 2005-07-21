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
  \brief Definitions of eigenvalue decomposition functions
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __eigen_h
#define __eigen_h

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>

namespace itpp {


  /*!
    \ingroup matrixdecomp
    \brief Calculates the eigenvalues and eigenvectors of a symmetric real matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real and symmetric \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]
    The eigenvectors are the columns of the matrix V.
    True is returned if the calculation was successful. Otherwise false.

    Uses the LAPACK routine DSYEV.
  */
  bool eig_sym(const mat &A, vec &d, mat &V);

  /*!
    \ingroup matrixdecomp
    \brief Calculates the eigenvalues of a symmetric real matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real and symmetric \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]
    True is returned if the calculation was successful. Otherwise false.

    Uses the LAPACK routine DSYEV.
  */
  bool eig_sym(const mat &A, vec &d);

  /*!
    \ingroup matrixdecomp
    \brief Calculates the eigenvalues of a symmetric real matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real and symmetric \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]

    Uses the LAPACK routine DSYEV.
  */
  vec eig_sym(const mat &A);

  /*!
    \ingroup matrixdecomp
    \brief Calculates the eigenvalues and eigenvectors of a hermitian complex matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex and hermitian \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]
    The eigenvectors are the columns of the matrix V.
    True is returned if the calculation was successful. Otherwise false.

    Uses the LAPACK routine ZHEEV.
  */
  bool eig_sym(const cmat &A, vec &d, cmat &V);

  /*!
    \ingroup matrixdecomp
    \brief Calculates the eigenvalues of a hermitian complex matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex and hermitian \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]
    True is returned if the calculation was successful. Otherwise false.

    Uses the LAPACK routine ZHEEV.
  */
  bool eig_sym(const cmat &A, vec &d);

  /*!
    \ingroup matrixdecomp
    \brief Calculates the eigenvalues of a hermitian complex matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex and hermitian \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]

    Uses the LAPACK routine ZHEEV.
  */
  vec eig_sym(const cmat &A);

  /*!
    \ingroup matrixdecomp
    \brief Caclulates the eigenvalues and eigenvectors of a real non-symmetric matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]
    The eigenvectors are the columns of the matrix V.
    True is returned if the calculation was successful. Otherwise false.

    Uses the LAPACK routine DGEEV.
  */
  bool eig(const mat &A, cvec &d, cmat &V);

  /*!
    \ingroup matrixdecomp
    \brief Caclulates the eigenvalues of a real non-symmetric matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]
    True is returned if the calculation was successful. Otherwise false.

    Uses the LAPACK routine DGEEV.
  */
  bool eig(const mat &A, cvec &d);

  /*!
    \ingroup matrixdecomp
    \brief Caclulates the eigenvalues of a real non-symmetric matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the real \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]

    Uses the LAPACK routine DGEEV.
  */
  cvec eig(const mat &A);

  /*!
    \ingroup matrixdecomp
    \brief Calculates the eigenvalues and eigenvectors of a complex non-hermitian matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]
    The eigenvectors are the columns of the matrix V.
    True is returned if the calculation was successful. Otherwise false.

    Uses the LAPACK routine ZGEEV.
  */
  bool eig(const cmat &A, cvec &d, cmat &V);

  /*!
    \ingroup matrixdecomp
    \brief Calculates the eigenvalues of a complex non-hermitian matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]
    True is returned if the calculation was successful. Otherwise false.

    Uses the LAPACK routine ZGEEV.
  */
  bool eig(const cmat &A, cvec &d);

  /*!
    \ingroup matrixdecomp
    \brief Calculates the eigenvalues of a complex non-hermitian matrix

    The Eigenvalues \f$\mathbf{d}(d_0, d_1, \ldots, d_{n-1})\f$ and the eigenvectors
    \f$\mathbf{v}_i, \: i=0, \ldots, n-1\f$ of the complex \f$n \times n\f$
    matrix \f$\mathbf{A}\f$ satisfies
    \f[
    \mathbf{A} \mathbf{v}_i = d_i \mathbf{v}_i\: i=0, \ldots, n-1.
    \f]

    Uses the LAPACK routine ZGEEV.
  */
  cvec eig(const cmat &A);


} //namespace itpp

#endif // __eigen_h
