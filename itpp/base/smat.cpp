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
  \brief Sparse Matrix Class implementation.
  \author Tony Ottosson and Tobias Ringström

  1.13

  2003/05/22 08:55:20
*/

#include <itpp/base/smat.h>

namespace itpp {
 

  /*--------------------------------------------------------------
   * Explicit initializations
   *--------------------------------------------------------------*/

  template class Sparse_Mat<int>;
  template class Sparse_Mat<double>;
  template class Sparse_Mat<std::complex<double> >;

  template sparse_imat operator+(const sparse_imat &, const sparse_imat &);
  template sparse_mat operator+(const sparse_mat &, const sparse_mat &);
  template sparse_cmat operator+(const sparse_cmat &, const sparse_cmat &);

  template sparse_imat operator*(const sparse_imat &, const sparse_imat &);
  template sparse_mat operator*(const sparse_mat &, const sparse_mat &);
  template sparse_cmat operator*(const sparse_cmat &, const sparse_cmat &);

  template ivec operator*(const ivec &, const sparse_imat &);
  template vec operator*(const vec &, const sparse_mat &);
  template cvec operator*(const cvec &, const sparse_cmat &);

  template ivec operator*(const sparse_imat &, const ivec &);
  template vec operator*(const sparse_mat &, const vec &);
  template cvec operator*(const sparse_cmat &, const cvec &);

  template imat trans_mult(const sparse_imat &);
  template mat trans_mult(const sparse_mat &);
  template cmat trans_mult(const sparse_cmat &);

  template sparse_imat trans_mult_s(const sparse_imat &);
  template sparse_mat trans_mult_s(const sparse_mat &);
  template sparse_cmat trans_mult_s(const sparse_cmat &);

  template sparse_imat trans_mult(const sparse_imat &, const sparse_imat &);
  template sparse_mat trans_mult(const sparse_mat &, const sparse_mat &);
  template sparse_cmat trans_mult(const sparse_cmat &, const sparse_cmat &);

  template ivec trans_mult(const sparse_imat &, const ivec &);
  template vec trans_mult(const sparse_mat &, const vec &);
  template cvec trans_mult(const sparse_cmat &, const cvec &);

  template sparse_imat mult_trans(const sparse_imat &, const sparse_imat &);
  template sparse_mat mult_trans(const sparse_mat &, const sparse_mat &);
  template sparse_cmat mult_trans(const sparse_cmat &, const sparse_cmat &);

} //namespace itpp
