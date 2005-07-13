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
  \brief Sparse Vector Class implementation
  \author Tony Ottosson and Tobias Ringström

  1.11

  2003/05/22 08:55:20
*/

#include "itpp/base/svec.h"

namespace itpp {

  /*--------------------------------------------------------------
   * Explicit initializations
   *--------------------------------------------------------------*/

  template class Sparse_Vec<int>;
  template class Sparse_Vec<double>;
  template class Sparse_Vec<std::complex<double> >;

  template sparse_ivec operator+(const sparse_ivec &, const sparse_ivec &);
  template sparse_vec operator+(const sparse_vec &, const sparse_vec &);
  template sparse_cvec operator+(const sparse_cvec &, const sparse_cvec &);

  template int operator*(const sparse_ivec &, const sparse_ivec &);
  template double operator*(const sparse_vec &, const sparse_vec &);
  template std::complex<double> operator*(const sparse_cvec &, const sparse_cvec &);

  template int operator*(const sparse_ivec &, const ivec &);
  template double operator*(const sparse_vec &, const vec &);
  template std::complex<double> operator*(const sparse_cvec &, const cvec &);

  template int operator*(const ivec &, const sparse_ivec &);
  template double operator*(const vec &, const sparse_vec &);
  template std::complex<double> operator*(const cvec &, const sparse_cvec &);

  template sparse_ivec elem_mult(const sparse_ivec &, const sparse_ivec &);
  template sparse_vec elem_mult(const sparse_vec &, const sparse_vec &);
  template sparse_cvec elem_mult(const sparse_cvec &, const sparse_cvec &);

  template ivec elem_mult(const sparse_ivec &, const ivec &);
  template vec elem_mult(const sparse_vec &, const vec &);
  template cvec elem_mult(const sparse_cvec &, const cvec &);

  template sparse_ivec elem_mult_s(const sparse_ivec &, const ivec &);
  template sparse_vec elem_mult_s(const sparse_vec &, const vec &);
  template sparse_cvec elem_mult_s(const sparse_cvec &, const cvec &);

  template ivec elem_mult(const ivec &, const sparse_ivec &);
  template vec elem_mult(const vec &, const sparse_vec &);
  template cvec elem_mult(const cvec &, const sparse_cvec &);

  template sparse_ivec elem_mult_s(const ivec &, const sparse_ivec &);
  template sparse_vec elem_mult_s(const vec &, const sparse_vec &);
  template sparse_cvec elem_mult_s(const cvec &, const sparse_cvec &);

} //namespace itpp
