/*!
 * \file
 * \brief Definitions of Singular Value Decompositions
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
 * Copyright (C) 1995-2006  (see AUTHORS file for a list of contributors)
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

#ifndef SVD_H
#define SVD_H

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


} // namespace itpp

#endif // #ifndef SVD_H
