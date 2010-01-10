/*!
 * \file
 * \brief Matrix Class Implementation
 * \author Tony Ottosson, Tobias Ringstrom, Adam Piatyszek and Conrad Sanderson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include <itpp/base/mat.h>

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#if defined (HAVE_BLAS)
#  include <itpp/base/blas.h>
#endif

//! \cond

namespace itpp
{

template<>
cmat cmat::hermitian_transpose() const
{
  cmat temp(no_cols, no_rows);
  for (int i = 0; i < no_rows; i++)
    for (int j = 0; j < no_cols; j++)
      temp(j, i) = std::conj(operator()(i,j));

  return temp;
}


// -------- Multiplication operator -------------

#if defined(HAVE_BLAS)

template<>
mat& mat::operator*=(const mat &m)
{
  it_assert_debug(no_cols == m.no_rows, "mat::operator*=(): Wrong sizes");
  mat r(no_rows, m.no_cols); // unnecessary memory??
  double alpha = 1.0;
  double beta = 0.0;
  char trans = 'n';
  blas::dgemm_(&trans, &trans, &no_rows, &m.no_cols, &no_cols, &alpha, data,
               &no_rows, m.data, &m.no_rows, &beta, r.data, &r.no_rows);
  operator=(r); // time consuming
  return *this;
}

template<>
cmat& cmat::operator*=(const cmat &m)
{
  it_assert_debug(no_cols == m.no_rows, "cmat::operator*=(): Wrong sizes");
  cmat r(no_rows, m.no_cols); // unnecessary memory??
  std::complex<double> alpha = std::complex<double>(1.0);
  std::complex<double> beta = std::complex<double>(0.0);
  char trans = 'n';
  blas::zgemm_(&trans, &trans, &no_rows, &m.no_cols, &no_cols, &alpha, data,
               &no_rows, m.data, &m.no_rows, &beta, r.data, &r.no_rows);
  operator=(r); // time consuming
  return *this;
}
#else
template<>
mat& mat::operator*=(const mat &m)
{
  it_assert_debug(no_cols == m.no_rows, "Mat<>::operator*=(): Wrong sizes");
  mat r(no_rows, m.no_cols);
  int r_pos = 0, pos = 0, m_pos = 0;

  for (int i = 0; i < r.no_cols; i++) {
    for (int j = 0; j < r.no_rows; j++) {
      double tmp = 0.0;
      pos = 0;
      for (int k = 0; k < no_cols; k++) {
        tmp += data[pos+j] * m.data[m_pos+k];
        pos += no_rows;
      }
      r.data[r_pos+j] = tmp;
    }
    r_pos += r.no_rows;
    m_pos += m.no_rows;
  }
  operator=(r); // time consuming
  return *this;
}

template<>
cmat& cmat::operator*=(const cmat &m)
{
  it_assert_debug(no_cols == m.no_rows, "Mat<>::operator*=(): Wrong sizes");
  cmat r(no_rows, m.no_cols);
  int r_pos = 0, pos = 0, m_pos = 0;

  for (int i = 0; i < r.no_cols; i++) {
    for (int j = 0; j < r.no_rows; j++) {
      std::complex<double> tmp(0.0);
      pos = 0;
      for (int k = 0; k < no_cols; k++) {
        tmp += data[pos+j] * m.data[m_pos+k];
        pos += no_rows;
      }
      r.data[r_pos+j] = tmp;
    }
    r_pos += r.no_rows;
    m_pos += m.no_rows;
  }
  operator=(r); // time consuming
  return *this;
}
#endif // HAVE_BLAS


#if defined(HAVE_BLAS)
template<>
mat operator*(const mat &m1, const mat &m2)
{
  it_assert_debug(m1.no_cols == m2.no_rows, "mat::operator*(): Wrong sizes");
  mat r(m1.no_rows, m2.no_cols);
  double alpha = 1.0;
  double beta = 0.0;
  char trans = 'n';
  blas::dgemm_(&trans, &trans, &m1.no_rows, &m2.no_cols, &m1.no_cols, &alpha,
               m1.data, &m1.no_rows, m2.data, &m2.no_rows, &beta, r.data,
               &r.no_rows);
  return r;
}

template<>
cmat operator*(const cmat &m1, const cmat &m2)
{
  it_assert_debug(m1.no_cols == m2.no_rows, "cmat::operator*(): Wrong sizes");
  cmat r(m1.no_rows, m2.no_cols);
  std::complex<double> alpha = std::complex<double>(1.0);
  std::complex<double> beta = std::complex<double>(0.0);
  char trans = 'n';
  blas::zgemm_(&trans, &trans, &m1.no_rows, &m2.no_cols, &m1.no_cols, &alpha,
               m1.data, &m1.no_rows, m2.data, &m2.no_rows, &beta, r.data,
               &r.no_rows);
  return r;
}
#else
template<>
mat operator*(const mat &m1, const mat &m2)
{
  it_assert_debug(m1.no_cols == m2.no_rows,
                  "Mat<>::operator*(): Wrong sizes");
  mat r(m1.no_rows, m2.no_cols);
  double *tr = r.data;
  double *t1;
  double *t2 = m2.data;
  for (int i = 0; i < r.no_cols; i++) {
    for (int j = 0; j < r.no_rows; j++) {
      double tmp = 0.0;
      t1 = m1.data + j;
      for (int k = m1.no_cols; k > 0; k--) {
        tmp += *(t1) * *(t2++);
        t1 += m1.no_rows;
      }
      *(tr++) = tmp;
      t2 -= m2.no_rows;
    }
    t2 += m2.no_rows;
  }
  return r;
}

template<>
cmat operator*(const cmat &m1, const cmat &m2)
{
  it_assert_debug(m1.no_cols == m2.no_rows,
                  "Mat<>::operator*(): Wrong sizes");
  cmat r(m1.no_rows, m2.no_cols);
  std::complex<double> *tr = r.data;
  std::complex<double> *t1;
  std::complex<double> *t2 = m2.data;
  for (int i = 0; i < r.no_cols; i++) {
    for (int j = 0; j < r.no_rows; j++) {
      std::complex<double> tmp(0.0);
      t1 = m1.data + j;
      for (int k = m1.no_cols; k > 0; k--) {
        tmp += *(t1) * *(t2++);
        t1 += m1.no_rows;
      }
      *(tr++) = tmp;
      t2 -= m2.no_rows;
    }
    t2 += m2.no_rows;
  }
  return r;
}
#endif // HAVE_BLAS


#if defined(HAVE_BLAS)
template<>
vec operator*(const mat &m, const vec &v)
{
  it_assert_debug(m.no_cols == v.size(), "mat::operator*(): Wrong sizes");
  vec r(m.no_rows);
  double alpha = 1.0;
  double beta = 0.0;
  char trans = 'n';
  int incr = 1;
  blas::dgemv_(&trans, &m.no_rows, &m.no_cols, &alpha, m.data, &m.no_rows,
               v._data(), &incr, &beta, r._data(), &incr);
  return r;
}

template<>
cvec operator*(const cmat &m, const cvec &v)
{
  it_assert_debug(m.no_cols == v.size(), "cmat::operator*(): Wrong sizes");
  cvec r(m.no_rows);
  std::complex<double> alpha = std::complex<double>(1.0);
  std::complex<double> beta = std::complex<double>(0.0);
  char trans = 'n';
  int incr = 1;
  blas::zgemv_(&trans, &m.no_rows, &m.no_cols, &alpha, m.data, &m.no_rows,
               v._data(), &incr, &beta, r._data(), &incr);
  return r;
}
#else
template<>
vec operator*(const mat &m, const vec &v)
{
  it_assert_debug(m.no_cols == v.size(),
                  "Mat<>::operator*(): Wrong sizes");
  vec r(m.no_rows);
  for (int i = 0; i < m.no_rows; i++) {
    r(i) = 0.0;
    int m_pos = 0;
    for (int k = 0; k < m.no_cols; k++) {
      r(i) += m.data[m_pos+i] * v(k);
      m_pos += m.no_rows;
    }
  }
  return r;
}

template<>
cvec operator*(const cmat &m, const cvec &v)
{
  it_assert_debug(m.no_cols == v.size(),
                  "Mat<>::operator*(): Wrong sizes");
  cvec r(m.no_rows);
  for (int i = 0; i < m.no_rows; i++) {
    r(i) = std::complex<double>(0.0);
    int m_pos = 0;
    for (int k = 0; k < m.no_cols; k++) {
      r(i) += m.data[m_pos+i] * v(k);
      m_pos += m.no_rows;
    }
  }
  return r;
}
#endif // HAVE_BLAS


//---------------------------------------------------------------------
// Instantiations
//---------------------------------------------------------------------

// class instantiations

template class Mat<double>;
template class Mat<std::complex<double> >;
template class Mat<int>;
template class Mat<short int>;
template class Mat<bin>;

// addition operators

template mat operator+(const mat &m1, const mat &m2);
template cmat operator+(const cmat &m1, const cmat &m2);
template imat operator+(const imat &m1, const imat &m2);
template smat operator+(const smat &m1, const smat &m2);
template bmat operator+(const bmat &m1, const bmat &m2);

template mat operator+(const mat &m, double t);
template cmat operator+(const cmat &m, std::complex<double> t);
template imat operator+(const imat &m, int t);
template smat operator+(const smat &m, short t);
template bmat operator+(const bmat &m, bin t);

template mat operator+(double t, const mat &m);
template cmat operator+(std::complex<double> t, const cmat &m);
template imat operator+(int t, const imat &m);
template smat operator+(short t, const smat &m);
template bmat operator+(bin t, const bmat &m);

// subraction operators

template mat operator-(const mat &m1, const mat &m2);
template cmat operator-(const cmat &m1, const cmat &m2);
template imat operator-(const imat &m1, const imat &m2);
template smat operator-(const smat &m1, const smat &m2);
template bmat operator-(const bmat &m1, const bmat &m2);

template mat operator-(const mat &m, double t);
template cmat operator-(const cmat &m, std::complex<double> t);
template imat operator-(const imat &m, int t);
template smat operator-(const smat &m, short t);
template bmat operator-(const bmat &m, bin t);

template mat operator-(double t, const mat &m);
template cmat operator-(std::complex<double> t, const cmat &m);
template imat operator-(int t, const imat &m);
template smat operator-(short t, const smat &m);
template bmat operator-(bin t, const bmat &m);

// unary minus

template mat operator-(const mat &m);
template cmat operator-(const cmat &m);
template imat operator-(const imat &m);
template smat operator-(const smat &m);
template bmat operator-(const bmat &m);

// multiplication operators

template imat operator*(const imat &m1, const imat &m2);
template smat operator*(const smat &m1, const smat &m2);
template bmat operator*(const bmat &m1, const bmat &m2);

template ivec operator*(const imat &m, const ivec &v);
template svec operator*(const smat &m, const svec &v);
template bvec operator*(const bmat &m, const bvec &v);

template mat operator*(const mat &m, double t);
template cmat operator*(const cmat &m, std::complex<double> t);
template imat operator*(const imat &m, int t);
template smat operator*(const smat &m, short t);
template bmat operator*(const bmat &m, bin t);

template mat operator*(double t, const mat &m);
template cmat operator*(std::complex<double> t, const cmat &m);
template imat operator*(int t, const imat &m);
template smat operator*(short t, const smat &m);
template bmat operator*(bin t, const bmat &m);

// elementwise multiplication

template mat elem_mult(const mat &m1, const mat &m2);
template cmat elem_mult(const cmat &m1, const cmat &m2);
template imat elem_mult(const imat &m1, const imat &m2);
template smat elem_mult(const smat &m1, const smat &m2);
template bmat elem_mult(const bmat &m1, const bmat &m2);

template void elem_mult_out(const mat &m1, const mat &m2, mat &out);
template void elem_mult_out(const cmat &m1, const cmat &m2, cmat &out);
template void elem_mult_out(const imat &m1, const imat &m2, imat &out);
template void elem_mult_out(const smat &m1, const smat &m2, smat &out);
template void elem_mult_out(const bmat &m1, const bmat &m2, bmat &out);

template void elem_mult_out(const mat &m1, const mat &m2,
                            const mat &m3, mat &out);
template void elem_mult_out(const cmat &m1, const cmat &m2,
                            const cmat &m3, cmat &out);
template void elem_mult_out(const imat &m1, const imat &m2,
                            const imat &m3, imat &out);
template void elem_mult_out(const smat &m1, const smat &m2,
                            const smat &m3, smat &out);
template void elem_mult_out(const bmat &m1, const bmat &m2,
                            const bmat &m3, bmat &out);

template void elem_mult_out(const mat &m1, const mat &m2, const mat &m3,
                            const mat &m4, mat &out);
template void elem_mult_out(const cmat &m1, const cmat &m2,
                            const cmat &m3, const cmat &m4, cmat &out);
template void elem_mult_out(const imat &m1, const imat &m2,
                            const imat &m3, const imat &m4, imat &out);
template void elem_mult_out(const smat &m1, const smat &m2,
                            const smat &m3, const smat &m4, smat &out);
template void elem_mult_out(const bmat &m1, const bmat &m2,
                            const bmat &m3, const bmat &m4, bmat &out);

template void elem_mult_inplace(const mat &m1, mat &m2);
template void elem_mult_inplace(const cmat &m1, cmat &m2);
template void elem_mult_inplace(const imat &m1, imat &m2);
template void elem_mult_inplace(const smat &m1, smat &m2);
template void elem_mult_inplace(const bmat &m1, bmat &m2);

template double elem_mult_sum(const mat &m1, const mat &m2);
template std::complex<double> elem_mult_sum(const cmat &m1, const cmat &m2);
template int elem_mult_sum(const imat &m1, const imat &m2);
template short elem_mult_sum(const smat &m1, const smat &m2);
template bin elem_mult_sum(const bmat &m1, const bmat &m2);

// division operator

template mat operator/(double t, const mat &m);
template cmat operator/(std::complex<double> t, const cmat &m);
template imat operator/(int t, const imat &m);
template smat operator/(short t, const smat &m);
template bmat operator/(bin t, const bmat &m);

template mat operator/(const mat &m, double t);
template cmat operator/(const cmat &m, std::complex<double> t);
template imat operator/(const imat &m, int t);
template smat operator/(const smat &m, short t);
template bmat operator/(const bmat &m, bin t);

// elementwise division

template mat elem_div(const mat &m1, const mat &m2);
template cmat elem_div(const cmat &m1, const cmat &m2);
template imat elem_div(const imat &m1, const imat &m2);
template smat elem_div(const smat &m1, const smat &m2);
template bmat elem_div(const bmat &m1, const bmat &m2);

template void elem_div_out(const mat &m1, const mat &m2, mat &out);
template void elem_div_out(const cmat &m1, const cmat &m2, cmat &out);
template void elem_div_out(const imat &m1, const imat &m2, imat &out);
template void elem_div_out(const smat &m1, const smat &m2, smat &out);
template void elem_div_out(const bmat &m1, const bmat &m2, bmat &out);

template double elem_div_sum(const mat &m1, const mat &m2);
template std::complex<double> elem_div_sum(const cmat &m1,
    const cmat &m2);
template int elem_div_sum(const imat &m1, const imat &m2);
template short elem_div_sum(const smat &m1, const smat &m2);
template bin elem_div_sum(const bmat &m1, const bmat &m2);

// concatenation

template mat concat_horizontal(const mat &m1, const mat &m2);
template cmat concat_horizontal(const cmat &m1, const cmat &m2);
template imat concat_horizontal(const imat &m1, const imat &m2);
template smat concat_horizontal(const smat &m1, const smat &m2);
template bmat concat_horizontal(const bmat &m1, const bmat &m2);

template mat concat_vertical(const mat &m1, const mat &m2);
template cmat concat_vertical(const cmat &m1, const cmat &m2);
template imat concat_vertical(const imat &m1, const imat &m2);
template smat concat_vertical(const smat &m1, const smat &m2);
template bmat concat_vertical(const bmat &m1, const bmat &m2);

// I/O streams

template std::ostream &operator<<(std::ostream &os, const mat  &m);
template std::ostream &operator<<(std::ostream &os, const cmat &m);
template std::ostream &operator<<(std::ostream &os, const imat  &m);
template std::ostream &operator<<(std::ostream &os, const smat  &m);
template std::ostream &operator<<(std::ostream &os, const bmat  &m);

template std::istream &operator>>(std::istream &is, mat  &m);
template std::istream &operator>>(std::istream &is, cmat &m);
template std::istream &operator>>(std::istream &is, imat  &m);
template std::istream &operator>>(std::istream &is, smat  &m);
template std::istream &operator>>(std::istream &is, bmat  &m);

} // namespace itpp

//! \endcond
