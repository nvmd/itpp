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
  \brief Definitions of QR factorisation functions
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __qr_h
#define __qr_h

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>

namespace itpp {


  /*! \addtogroup matrixdecomp
   */
  //!@{
  /*!
    \brief QR factorisation of real matrix

    The QR factorization of the real matrix \f$\mathbf{A}\f$ of size \f$m \times n\f$ is given
    by
    \f[
    \mathbf{A} = \mathbf{Q} \mathbf{R} ,
    \f]
    where \f$\mathbf{Q}\f$ is an \f$m \times m\f$ orthogonal matrix and \f$\mathbf{R}\f$ is an \f$m \times n\f$ upper triangular matrix.

    Returns true is calculation succeeds. False otherwise.
    Uses the LAPACK routine DGEQRF and DORGQR.
  */
  bool qr(const mat &A, mat &Q, mat &R);

  /*!
    \brief QR factorisation of real matrix with pivoting

    The QR factorization of the real matrix \f$\mathbf{A}\f$ of size \f$m \times n\f$ is given
    by
    \f[
    \mathbf{A} \mathbf{P} = \mathbf{Q} \mathbf{R} ,
    \f]
    where \f$\mathbf{Q}\f$ is an \f$m \times m\f$ orthogonal matrix, \f$\mathbf{R}\f$ is an \f$m \times n\f$ upper triangular matrix
    and \f$\mathbf{P}\f$ is an \f$n \times n\f$ permutation matrix.

    Returns true is calculation succeeds. False otherwise.
    Uses the LAPACK routines DGEQP3 and DORGQR.
  */
  bool qr(const mat &A, mat &Q, mat &R, bmat &P);

  /*!
    \brief QR factorisation of a complex matrix

    The QR factorization of the complex matrix \f$\mathbf{A}\f$ of size \f$m \times n\f$ is given
    by
    \f[
    \mathbf{A} = \mathbf{Q} \mathbf{R} ,
    \f]
    where \f$\mathbf{Q}\f$ is an \f$m \times m\f$ unitary matrix and \f$\mathbf{R}\f$ is an \f$m \times n\f$ upper triangular matrix.

    Returns true is calculation succeeds. False otherwise.
    Uses the LAPACK routines ZGEQRF and ZUNGQR.
  */
  bool qr(const cmat &A, cmat &Q, cmat &R);

  /*!
    \brief QR factorisation of a complex matrix with pivoting

    The QR factorization of the complex matrix \f$\mathbf{A}\f$ of size \f$m \times n\f$ is given
    by
    \f[
    \mathbf{A} \mathbf{P} = \mathbf{Q} \mathbf{R} ,
    \f]
    where \f$\mathbf{Q}\f$ is an \f$m \times m\f$ unitary matrix, \f$\mathbf{R}\f$ is an \f$m \times n\f$ upper triangular matrix
    and \f$\mathbf{P}\f$ is an \f$n \times n\f$ permutation matrix.

    Returns true is calculation succeeds. False otherwise.
    Uses the LAPACK routines ZGEQP3 and ZUNGQR.
  */
  bool qr(const cmat &A, cmat &Q, cmat &R, bmat &P);

  //!@}


} //namespace itpp

#endif // __qr_h
