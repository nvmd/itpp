/*!
 * \file
 * \brief Sparse Matrix Class implementation
 * \author Tony Ottosson and Tobias Ringstrom
 * 
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
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
 * -------------------------------------------------------------------------
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

} // namespace itpp
