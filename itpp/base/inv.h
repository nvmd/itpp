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
  \brief Inverse of square matrices
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __inv_h
#define __inv_h

#include <itpp/base/mat.h>

namespace itpp {

  /*!
    \addtogroup inverse
  */

  /*!
    \brief Inverse of real square matrix.
    \ingroup inverse

    Calculate the inverse of the real matrix \f$\mathbf{X}\f$

    Solves the equation system \f$ \mathbf{Y} \mathbf{X} = \mathbf{I}\f$ using LU-factorization.
    IT++ needs to be compiled with the LAPACK for the inverse to be available.
  */
  bool inv(const mat &X, mat &Y);

  /*!
    \brief Inverse of real square matrix.
    \ingroup inverse

    Calculate the inverse of the real matrix \f$\mathbf{X}\f$

    Solves the equation system \f$ \mathbf{Y} \mathbf{X} = \mathbf{I}\f$ using LU-factorization.
    IT++ needs to be compiled with LAPACK support for the inverse to be available.
  */
  mat inv(const mat &X);


  /*!
    \brief Inverse of complex square matrix.
    \ingroup inverse

    Calculate the inverse of the complex matrix \f$\mathbf{X}\f$

    Solves the equation system \f$ \mathbf{Y} \mathbf{X} = \mathbf{I}\f$ using LU-factorization.
    IT++ needs to be compiled with LAPACK support for the inverse to be available.
  */
  bool inv(const cmat &X, cmat &Y);

  /*!
    \brief Inverse of real square matrix.
    \ingroup inverse

    Calculate the inverse of the complex matrix \f$\mathbf{X}\f$

    Solves the equation system \f$ \mathbf{Y} \mathbf{X} = \mathbf{I}\f$ using LU-factorization.
    IT++ needs to be compiled with the LAPACK for the inverse to be available.
  */
  cmat inv(const cmat &X);


} //namespace itpp

#endif // __inv_h
