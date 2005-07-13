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
  \brief Definitions of Cholesky factorisation functions
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __cholesky_h
#define __cholesky_h

#include "itpp/base/vec.h"
#include "itpp/base/mat.h"

namespace itpp {


  /*! \addtogroup matrixdecomp
   */
  //!@{

  /*!
    \brief Cholesky factorisation of real symmetric and positive definite matrix

    The Cholesky factorisation of a real symmetric positive-definite matrix \f$\mathbf{X}\f$
    of size \f$n \times n\f$ is given by
    \f[
    \mathbf{X} = \mathbf{F}^T \mathbf{F}
    \f]
    where \f$\mathbf{F}\f$ is an upper trangular \f$n \times n\f$ matrix.

    Returns true if calcuation succeeded. False otherwise.
  */
  bool chol(const mat &X, mat &F);

  /*!
    \brief Cholesky factorisation of real symmetric and positive definite matrix

    The Cholesky factorisation of a real symmetric positive-definite matrix \f$\mathbf{X}\f$
    of size \f$n \times n\f$ is given by
    \f[
    \mathbf{X} = \mathbf{F}^T \mathbf{F}
    \f]
    where \f$\mathbf{F}\f$ is an upper trangular \f$n \times n\f$ matrix.
  */
  mat chol(const mat &X);


  /*!
    \brief Cholesky factorisation of complex hermitian and positive-definite matrix

    The Cholesky factorisation of a hermitian positive-definite matrix \f$\mathbf{X}\f$
    of size \f$n \times n\f$ is given by
    \f[
    \mathbf{X} = \mathbf{F}^H \mathbf{F}
    \f]
    where \f$\mathbf{F}\f$ is an upper trangular \f$n \times n\f$ matrix.

    Returns true if calcuation succeeded. False otherwise.

    If \c X is positive definite, true is returned and \c F=chol(X)
    produces an upper triangular \c F. If also \c X is symmetric then \c F'*F = X.
    If \c X is not positive definite, false is returned.
  */
  bool chol(const cmat &X, cmat &F);

  /*!
    \brief Cholesky factorisation of complex hermitian and positive-definite matrix

    The Cholesky factorisation of a hermitian positive-definite matrix \f$\mathbf{X}\f$
    of size \f$n \times n\f$ is given by
    \f[
    \mathbf{X} = \mathbf{F}^H \mathbf{F}
    \f]
    where \f$\mathbf{F}\f$ is an upper trangular \f$n \times n\f$ matrix.
  */
  cmat chol(const cmat &X);

  //!@}

} //namespace itpp

#endif // __cholesky_h
