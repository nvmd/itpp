/*!
 * \file
 * \brief Definitions of QR factorisation functions
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
 * Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

#ifndef QR_H
#define QR_H

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


} // namespace itpp

#endif // #ifndef QR_H
