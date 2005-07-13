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
  \brief Definitions of LU factorisation functions
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __lu_h
#define __lu_h

#include "itpp/base/vec.h"
#include "itpp/base/mat.h"

namespace itpp {


  /*! \addtogroup matrixdecomp
   */
  //!@{
  /*!
    \brief LU factorisation of real matrix

    The LU factorization of the real matrix \f$\mathbf{X}\f$ of size \f$n \times n\f$ is given
    by
    \f[
    \mathbf{X} = \mathbf{P}^T \mathbf{L} \mathbf{U} ,
    \f]
    where \f$\mathbf{L}\f$ and \f$\mathbf{U}\f$ are lower and upper triangular matrices
    and \f$\mathbf{P}\f$ is a permutation matrix.

    The interchange permutation vector \a p is such that \a k and \a p(k) should be
    changed for all \a k. Given this vector a permuation matrix can be constructed using the
    function
    \code
    bmat permuation_matrix(const ivec &p)
    \endcode

    If \a X is an \a n by \a n matrix \a lu(X,L,U,p) computes the LU decomposition.
    \a L is a lower trangular, \a U an upper triangular matrix.
    \a p is the interchange permutation vector such that \a k and \a p(k) should be
    changed for all \a k.

    Returns true is calculation succeeds. False otherwise.
  */
  bool lu(const mat &X, mat &L, mat &U, ivec &p);


  /*!
    \brief LU factorisation of real matrix

    The LU factorization of the complex matrix \f$\mathbf{X}\f$ of size \f$n \times n\f$ is given
    by
    \f[
    \mathbf{X} = \mathbf{P}^T \mathbf{L} \mathbf{U} ,
    \f]
    where \f$\mathbf{L}\f$ and \f$\mathbf{U}\f$ are lower and upper triangular matrices
    and \f$\mathbf{P}\f$ is a permutation matrix.

    The interchange permutation vector \a p is such that \a k and \a p(k) should be
    changed for all \a k. Given this vector a permuation matrix can be constructed using the
    function
    \code
    bmat permuation_matrix(const ivec &p)
    \endcode

    If \a X is an \a n by \a n matrix \a lu(X,L,U,p) computes the LU decomposition.
    \a L is a lower trangular, \a U an upper triangular matrix.
    \a p is the interchange permutation vector such that elements \a k and row \a p(k) should be
    interchanged.

    Returns true is calculation succeeds. False otherwise.
  */
  bool lu(const cmat &X, cmat &L, cmat &U, ivec &p);


  //! Makes swapping of vector b according to the inerchange permutation vector p.
  void interchange_permutations(vec &b, const ivec &p);

  //! Make permutation matrix P from the interchange permutation vector p.
  bmat permutation_matrix(const ivec &p);
  //!@}

} //namespace itpp

#endif // __lu_h
