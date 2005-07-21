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
  \brief Definitions of Singular Value Decompositions
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __svd_h
#define __svd_h

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>

namespace itpp {

  /*!
    \ingroup matrixdecomp
    \brief Singular Value Decomposition (SVD)

    The svd-algorithm computes the decomposition of a real \f$m \times n\f$ matrix \f$\mathbf{A}\f$
    so that
    \f[
    \mathbf{U}^T \mathbf{A} \mathbf{V} = \mathrm{diag}(\mathbf{s}) = \mathrm{diag}(\sigma_1, \ldots, \sigma_p)
    \f]
    where the elements of \f$\mathbf{s}\f$, \f$\sigma_1 \geq \sigma_2 \geq \ldots \sigma_p \geq 0\f$
    are the singular values of \f$\mathbf{A}\f$. Or put differently
    \f[
    \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^H
    \f]
  */
  bool svd(const mat &A, vec &S);

  /*!
    \ingroup matrixdecomp
    \brief Singular Value Decomposition (SVD)

    The svd-algorithm computes the decomposition of a complex \f$m \times n\f$ matrix \f$\mathbf{A}\f$
    so that
    \f[
    \mathbf{U}^T \mathbf{A} \mathbf{V} = \mathrm{diag}(\mathbf{s}) = \mathrm{diag}(\sigma_1, \ldots, \sigma_p)
    \f]
    where the elements of \f$\mathbf{s}\f$, \f$\sigma_1 \geq \sigma_2 \geq \ldots \sigma_p \geq 0\f$
    are the singular values of \f$\mathbf{A}\f$. Or put differently
    \f[
    \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^H
    \f]
  */
  bool svd(const cmat &A, vec &S);

  /*!
    \ingroup matrixdecomp
    \brief Singular Value Decomposition (SVD)

    The svd-algorithm computes the decomposition of a real \f$m \times n\f$ matrix \f$\mathbf{A}\f$
    so that
    \f[
    \mathbf{U}^T \mathbf{A} \mathbf{V} = \mathrm{diag}(\mathbf{s}) = \mathrm{diag}(\sigma_1, \ldots, \sigma_p)
    \f]
    where the elements of \f$\mathbf{s}\f$, \f$\sigma_1 \geq \sigma_2 \geq \ldots \sigma_p \geq 0\f$
    are the singular values of \f$\mathbf{A}\f$. Or put differently
    \f[
    \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^H
    \f]
  */
  vec svd(const mat &A);

  /*!
    \ingroup matrixdecomp
    \brief Singular Value Decomposition (SVD)

    The svd-algorithm computes the decomposition of a complex \f$m \times n\f$ matrix \f$\mathbf{A}\f$
    so that
    \f[
    \mathbf{U}^T \mathbf{A} \mathbf{V} = \mathrm{diag}(\mathbf{s}) = \mathrm{diag}(\sigma_1, \ldots, \sigma_p)
    \f]
    where the elements of \f$\mathbf{s}\f$, \f$\sigma_1 \geq \sigma_2 \geq \ldots \sigma_p \geq 0\f$
    are the singular values of \f$\mathbf{A}\f$. Or put differently
    \f[
    \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^H
    \f]
  */
  vec svd(const cmat &A);

  /*!
    \ingroup matrixdecomp
    \brief Singular Value Decomposition (SVD)

    The svd-algorithm computes the decomposition of a real \f$m \times n\f$ matrix \f$\mathbf{A}\f$
    so that
    \f[
    \mathbf{U}^T \mathbf{A} \mathbf{V} = \mathrm{diag}(\mathbf{s}) = \mathrm{diag}(\sigma_1, \ldots, \sigma_p)
    \f]
    where the elements of \f$\mathbf{s}\f$, \f$\sigma_1 \geq \sigma_2 \geq \ldots \sigma_p \geq 0\f$
    are the singular values of \f$\mathbf{A}\f$. Or put differently
    \f[
    \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^H
    \f]
  */
  bool svd(const mat &A, mat &U, vec &S, mat &V);

  /*!
    \ingroup matrixdecomp
    \brief Singular Value Decomposition (SVD)

    The svd-algorithm computes the decomposition of a complex \f$m \times n\f$ matrix \f$\mathbf{A}\f$
    so that
    \f[
    \mathbf{U}^T \mathbf{A} \mathbf{V} = \mathrm{diag}(\mathbf{s}) = \mathrm{diag}(\sigma_1, \ldots, \sigma_p)
    \f]
    where the elements of \f$\mathbf{s}\f$, \f$\sigma_1 \geq \sigma_2 \geq \ldots \sigma_p \geq 0\f$
    are the singular values of \f$\mathbf{A}\f$. Or put differently
    \f[
    \mathbf{A} = \mathbf{U} \mathbf{S} \mathbf{V}^H
    \f]
  */
  bool svd(const cmat &A, cmat &U, vec &S, cmat &V);


} //namespace itpp

#endif // __svd_h
