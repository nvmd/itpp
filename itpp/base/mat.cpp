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

template<>
mat& mat::operator*=(const mat &m)
{
  it_assert_debug(no_cols == m.no_rows, "mat::operator*=(): Wrong sizes");
  mat r(no_rows, m.no_cols); // unnecessary memory??
#if defined(HAVE_BLAS)
  double alpha = 1.0;
  double beta = 0.0;
  char trans = 'n';
  int_blas_t _no_rows   = static_cast<int_blas_t>(no_rows);
  int_blas_t _no_rows_m = static_cast<int_blas_t>(m.no_rows);
  int_blas_t _no_rows_r = static_cast<int_blas_t>(r.no_rows);
  int_blas_t _no_cols   = static_cast<int_blas_t>(no_cols);
  int_blas_t _no_cols_m = static_cast<int_blas_t>(m.no_cols);
  blas::dgemm_(&trans, &trans, &_no_rows, &_no_cols_m, &_no_cols, &alpha, data,
               &_no_rows, m.data, &_no_rows_m, &beta, r.data, &_no_rows_r);
#else
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
#endif
  operator=(r); // time consuming
  return *this;
}

template<>
cmat& cmat::operator*=(const cmat &m)
{
  it_assert_debug(no_cols == m.no_rows, "cmat::operator*=(): Wrong sizes");
  cmat r(no_rows, m.no_cols); // unnecessary memory??
#if defined(HAVE_BLAS)
  std::complex<double> alpha = std::complex<double>(1.0);
  std::complex<double> beta = std::complex<double>(0.0);
  char trans = 'n';
  int_blas_t _no_rows   = static_cast<int_blas_t>(no_rows);
  int_blas_t _no_rows_m = static_cast<int_blas_t>(m.no_rows);
  int_blas_t _no_rows_r = static_cast<int_blas_t>(r.no_rows);
  int_blas_t _no_cols   = static_cast<int_blas_t>(no_cols);
  int_blas_t _no_cols_m = static_cast<int_blas_t>(m.no_cols);
  blas::zgemm_(&trans, &trans, &_no_rows, &_no_cols_m, &_no_cols, &alpha, data,
               &_no_rows, m.data, &_no_rows_m, &beta, r.data, &_no_rows_r);
#else
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
#endif
  operator=(r); // time consuming
  return *this;
}

template<>
mat operator*(const mat &m1, const mat &m2)
{
  it_assert_debug(m1.cols() == m2.rows(), "mat::operator*(): Wrong sizes");
  mat r(m1.rows(), m2.cols());
#if defined(HAVE_BLAS)
  double alpha = 1.0;
  double beta = 0.0;
  char trans = 'n';
  int_blas_t m1_r = static_cast<int_blas_t>(m1.rows());
  int_blas_t m1_c = static_cast<int_blas_t>(m1.cols());
  int_blas_t m2_r = static_cast<int_blas_t>(m2.rows());
  int_blas_t m2_c = static_cast<int_blas_t>(m2.cols());
  blas::dgemm_(&trans, &trans, &m1_r, &m2_c, &m1_c, &alpha,
               m1._data(), &m1_r, m2._data(), &m2_r, &beta, r._data(),
               &m1_r);
#else
  double *tr = r._data();
  const double *t1;
  const double *t2 = m2._data();
  for (int i = 0; i < r.cols(); i++) {
    for (int j = 0; j < r.rows(); j++) {
      double tmp = 0.0;
      t1 = m1._data() + j;
      for (int k = m1.cols(); k > 0; k--) {
        tmp += *(t1) * *(t2++);
        t1 += m1.rows();
      }
      *(tr++) = tmp;
      t2 -= m2.rows();
    }
    t2 += m2.rows();
  }
#endif
  return r;
}

template<>
cmat operator*(const cmat &m1, const cmat &m2)
{
  it_assert_debug(m1.cols() == m2.rows(), "cmat::operator*(): Wrong sizes");
  cmat r(m1.rows(), m2.cols());
#if defined(HAVE_BLAS)
  std::complex<double> alpha = std::complex<double>(1.0);
  std::complex<double> beta = std::complex<double>(0.0);
  char trans = 'n';
  int_blas_t m1_r = static_cast<int_blas_t>(m1.rows());
  int_blas_t m1_c = static_cast<int_blas_t>(m1.cols());
  int_blas_t m2_r = static_cast<int_blas_t>(m2.rows());
  int_blas_t m2_c = static_cast<int_blas_t>(m2.cols());
  blas::zgemm_(&trans, &trans, &m1_r, &m2_c, &m1_c, &alpha,
               m1._data(), &m1_r, m2._data(), &m2_r, &beta, r._data(),
               &m1_r);
#else
  std::complex<double> *tr = r._data();
  const std::complex<double> *t1;
  const std::complex<double> *t2 = m2._data();
  for (int i = 0; i < r.cols(); i++) {
    for (int j = 0; j < r.rows(); j++) {
      std::complex<double> tmp(0.0);
      t1 = m1._data() + j;
      for (int k = m1.cols(); k > 0; k--) {
        tmp += *(t1) * *(t2++);
        t1 += m1.rows();
      }
      *(tr++) = tmp;
      t2 -= m2.rows();
    }
    t2 += m2.rows();
  }
#endif
  return r;
}

template<>
vec operator*(const mat &m, const vec &v)
{
  it_assert_debug(m.cols() == v.size(), "mat::operator*(): Wrong sizes");
  vec r(m.rows());
#if defined(HAVE_BLAS)
  double alpha = 1.0;
  double beta = 0.0;
  char trans = 'n';
  int_blas_t incr = 1;
  int_blas_t m_r = static_cast<int_blas_t>(m.rows());
  int_blas_t m_c = static_cast<int_blas_t>(m.cols());
  blas::dgemv_(&trans, &m_r, &m_c, &alpha, m._data(), &m_r,
               v._data(), &incr, &beta, r._data(), &incr);
#else
  for (int i = 0; i < m.rows(); i++) {
    r(i) = 0.0;
    int m_pos = 0;
    for (int k = 0; k < m.cols(); k++) {
      r(i) += m._data()[m_pos+i] * v(k);
      m_pos += m.rows();
    }
  }
#endif
  return r;
}

template<>
cvec operator*(const cmat &m, const cvec &v)
{
  it_assert_debug(m.cols() == v.size(), "cmat::operator*(): Wrong sizes");
  cvec r(m.rows());
#if defined(HAVE_BLAS)
  std::complex<double> alpha = std::complex<double>(1.0);
  std::complex<double> beta = std::complex<double>(0.0);
  char trans = 'n';
  int_blas_t incr = 1;
  int_blas_t m_r = static_cast<int_blas_t>(m.rows());
  int_blas_t m_c = static_cast<int_blas_t>(m.cols());
  blas::zgemv_(&trans, &m_r, &m_c, &alpha, m._data(), &m_r,
               v._data(), &incr, &beta, r._data(), &incr);
#else
  for (int i = 0; i < m.rows(); i++) {
    r(i) = std::complex<double>(0.0);
    int m_pos = 0;
    for (int k = 0; k < m.cols(); k++) {
      r(i) += m._data()[m_pos+i] * v(k);
      m_pos += m.rows();
    }
  }
#endif
  return r;
}


//---------------------------------------------------------------------
// Instantiations
//---------------------------------------------------------------------

// class instantiations

template class ITPP_EXPORT Mat<double>;
template class ITPP_EXPORT Mat<std::complex<double> >;
template class ITPP_EXPORT Mat<int>;
template class ITPP_EXPORT Mat<short int>;
template class ITPP_EXPORT Mat<bin>;

// addition operators

template ITPP_EXPORT  mat operator+(const mat &m1, const mat &m2);
template ITPP_EXPORT  cmat operator+(const cmat &m1, const cmat &m2);
template ITPP_EXPORT  imat operator+(const imat &m1, const imat &m2);
template ITPP_EXPORT  smat operator+(const smat &m1, const smat &m2);
template ITPP_EXPORT  bmat operator+(const bmat &m1, const bmat &m2);

template ITPP_EXPORT  mat operator+(const mat &m, double t);
template ITPP_EXPORT  cmat operator+(const cmat &m, std::complex<double> t);
template ITPP_EXPORT  imat operator+(const imat &m, int t);
template ITPP_EXPORT  smat operator+(const smat &m, short t);
template ITPP_EXPORT  bmat operator+(const bmat &m, bin t);

template ITPP_EXPORT  mat operator+(double t, const mat &m);
template ITPP_EXPORT  cmat operator+(std::complex<double> t, const cmat &m);
template ITPP_EXPORT  imat operator+(int t, const imat &m);
template ITPP_EXPORT  smat operator+(short t, const smat &m);
template ITPP_EXPORT  bmat operator+(bin t, const bmat &m);

// subraction operators

template ITPP_EXPORT  mat operator-(const mat &m1, const mat &m2);
template ITPP_EXPORT  cmat operator-(const cmat &m1, const cmat &m2);
template ITPP_EXPORT  imat operator-(const imat &m1, const imat &m2);
template ITPP_EXPORT  smat operator-(const smat &m1, const smat &m2);
template ITPP_EXPORT  bmat operator-(const bmat &m1, const bmat &m2);

template ITPP_EXPORT  mat operator-(const mat &m, double t);
template ITPP_EXPORT  cmat operator-(const cmat &m, std::complex<double> t);
template ITPP_EXPORT  imat operator-(const imat &m, int t);
template ITPP_EXPORT  smat operator-(const smat &m, short t);
template ITPP_EXPORT  bmat operator-(const bmat &m, bin t);

template ITPP_EXPORT  mat operator-(double t, const mat &m);
template ITPP_EXPORT  cmat operator-(std::complex<double> t, const cmat &m);
template ITPP_EXPORT  imat operator-(int t, const imat &m);
template ITPP_EXPORT  smat operator-(short t, const smat &m);
template ITPP_EXPORT  bmat operator-(bin t, const bmat &m);

// unary minus

template ITPP_EXPORT  mat operator-(const mat &m);
template ITPP_EXPORT  cmat operator-(const cmat &m);
template ITPP_EXPORT  imat operator-(const imat &m);
template ITPP_EXPORT  smat operator-(const smat &m);
template ITPP_EXPORT  bmat operator-(const bmat &m);

// multiplication operators

template ITPP_EXPORT  imat operator*(const imat &m1, const imat &m2);
template ITPP_EXPORT  smat operator*(const smat &m1, const smat &m2);
template ITPP_EXPORT  bmat operator*(const bmat &m1, const bmat &m2);

template ITPP_EXPORT  ivec operator*(const imat &m, const ivec &v);
template ITPP_EXPORT  svec operator*(const smat &m, const svec &v);
template ITPP_EXPORT  bvec operator*(const bmat &m, const bvec &v);

template ITPP_EXPORT  mat operator*(const mat &m, double t);
template ITPP_EXPORT  cmat operator*(const cmat &m, std::complex<double> t);
template ITPP_EXPORT  imat operator*(const imat &m, int t);
template ITPP_EXPORT  smat operator*(const smat &m, short t);
template ITPP_EXPORT  bmat operator*(const bmat &m, bin t);

template ITPP_EXPORT  mat operator*(double t, const mat &m);
template ITPP_EXPORT  cmat operator*(std::complex<double> t, const cmat &m);
template ITPP_EXPORT  imat operator*(int t, const imat &m);
template ITPP_EXPORT  smat operator*(short t, const smat &m);
template ITPP_EXPORT  bmat operator*(bin t, const bmat &m);

// elementwise multiplication

template ITPP_EXPORT  mat elem_mult(const mat &m1, const mat &m2);
template ITPP_EXPORT  cmat elem_mult(const cmat &m1, const cmat &m2);
template ITPP_EXPORT  imat elem_mult(const imat &m1, const imat &m2);
template ITPP_EXPORT  smat elem_mult(const smat &m1, const smat &m2);
template ITPP_EXPORT  bmat elem_mult(const bmat &m1, const bmat &m2);

template ITPP_EXPORT  void elem_mult_out(const mat &m1, const mat &m2, mat &out);
template ITPP_EXPORT  void elem_mult_out(const cmat &m1, const cmat &m2, cmat &out);
template ITPP_EXPORT  void elem_mult_out(const imat &m1, const imat &m2, imat &out);
template ITPP_EXPORT  void elem_mult_out(const smat &m1, const smat &m2, smat &out);
template ITPP_EXPORT  void elem_mult_out(const bmat &m1, const bmat &m2, bmat &out);

template ITPP_EXPORT  void elem_mult_out(const mat &m1, const mat &m2,
                            const mat &m3, mat &out);
template ITPP_EXPORT  void elem_mult_out(const cmat &m1, const cmat &m2,
                            const cmat &m3, cmat &out);
template ITPP_EXPORT  void elem_mult_out(const imat &m1, const imat &m2,
                            const imat &m3, imat &out);
template ITPP_EXPORT  void elem_mult_out(const smat &m1, const smat &m2,
                            const smat &m3, smat &out);
template ITPP_EXPORT  void elem_mult_out(const bmat &m1, const bmat &m2,
                            const bmat &m3, bmat &out);

template ITPP_EXPORT  void elem_mult_out(const mat &m1, const mat &m2, const mat &m3,
                            const mat &m4, mat &out);
template ITPP_EXPORT  void elem_mult_out(const cmat &m1, const cmat &m2,
                            const cmat &m3, const cmat &m4, cmat &out);
template ITPP_EXPORT  void elem_mult_out(const imat &m1, const imat &m2,
                            const imat &m3, const imat &m4, imat &out);
template ITPP_EXPORT  void elem_mult_out(const smat &m1, const smat &m2,
                            const smat &m3, const smat &m4, smat &out);
template ITPP_EXPORT  void elem_mult_out(const bmat &m1, const bmat &m2,
                            const bmat &m3, const bmat &m4, bmat &out);

template ITPP_EXPORT  void elem_mult_inplace(const mat &m1, mat &m2);
template ITPP_EXPORT  void elem_mult_inplace(const cmat &m1, cmat &m2);
template ITPP_EXPORT  void elem_mult_inplace(const imat &m1, imat &m2);
template ITPP_EXPORT  void elem_mult_inplace(const smat &m1, smat &m2);
template ITPP_EXPORT  void elem_mult_inplace(const bmat &m1, bmat &m2);

template ITPP_EXPORT  double elem_mult_sum(const mat &m1, const mat &m2);
template ITPP_EXPORT  std::complex<double> elem_mult_sum(const cmat &m1, const cmat &m2);
template ITPP_EXPORT  int elem_mult_sum(const imat &m1, const imat &m2);
template ITPP_EXPORT  short elem_mult_sum(const smat &m1, const smat &m2);
template ITPP_EXPORT  bin elem_mult_sum(const bmat &m1, const bmat &m2);

// division operator

template ITPP_EXPORT  mat operator/(double t, const mat &m);
template ITPP_EXPORT  cmat operator/(std::complex<double> t, const cmat &m);
template ITPP_EXPORT  imat operator/(int t, const imat &m);
template ITPP_EXPORT  smat operator/(short t, const smat &m);
template ITPP_EXPORT  bmat operator/(bin t, const bmat &m);

template ITPP_EXPORT  mat operator/(const mat &m, double t);
template ITPP_EXPORT  cmat operator/(const cmat &m, std::complex<double> t);
template ITPP_EXPORT  imat operator/(const imat &m, int t);
template ITPP_EXPORT  smat operator/(const smat &m, short t);
template ITPP_EXPORT  bmat operator/(const bmat &m, bin t);

// elementwise division

template ITPP_EXPORT  mat elem_div(const mat &m1, const mat &m2);
template ITPP_EXPORT  cmat elem_div(const cmat &m1, const cmat &m2);
template ITPP_EXPORT  imat elem_div(const imat &m1, const imat &m2);
template ITPP_EXPORT  smat elem_div(const smat &m1, const smat &m2);
template ITPP_EXPORT  bmat elem_div(const bmat &m1, const bmat &m2);

template ITPP_EXPORT  void elem_div_out(const mat &m1, const mat &m2, mat &out);
template ITPP_EXPORT  void elem_div_out(const cmat &m1, const cmat &m2, cmat &out);
template ITPP_EXPORT  void elem_div_out(const imat &m1, const imat &m2, imat &out);
template ITPP_EXPORT  void elem_div_out(const smat &m1, const smat &m2, smat &out);
template ITPP_EXPORT  void elem_div_out(const bmat &m1, const bmat &m2, bmat &out);

template ITPP_EXPORT  double elem_div_sum(const mat &m1, const mat &m2);
template ITPP_EXPORT  std::complex<double> elem_div_sum(const cmat &m1,
    const cmat &m2);
template ITPP_EXPORT  int elem_div_sum(const imat &m1, const imat &m2);
template ITPP_EXPORT  short elem_div_sum(const smat &m1, const smat &m2);
template ITPP_EXPORT  bin elem_div_sum(const bmat &m1, const bmat &m2);

// concatenation

template ITPP_EXPORT  mat concat_horizontal(const mat &m1, const mat &m2);
template ITPP_EXPORT  cmat concat_horizontal(const cmat &m1, const cmat &m2);
template ITPP_EXPORT  imat concat_horizontal(const imat &m1, const imat &m2);
template ITPP_EXPORT  smat concat_horizontal(const smat &m1, const smat &m2);
template ITPP_EXPORT  bmat concat_horizontal(const bmat &m1, const bmat &m2);

template ITPP_EXPORT  mat concat_vertical(const mat &m1, const mat &m2);
template ITPP_EXPORT  cmat concat_vertical(const cmat &m1, const cmat &m2);
template ITPP_EXPORT  imat concat_vertical(const imat &m1, const imat &m2);
template ITPP_EXPORT  smat concat_vertical(const smat &m1, const smat &m2);
template ITPP_EXPORT  bmat concat_vertical(const bmat &m1, const bmat &m2);

// I/O streams

template ITPP_EXPORT  std::ostream &operator<<(std::ostream &os, const mat  &m);
template ITPP_EXPORT  std::ostream &operator<<(std::ostream &os, const cmat &m);
template ITPP_EXPORT  std::ostream &operator<<(std::ostream &os, const imat  &m);
template ITPP_EXPORT  std::ostream &operator<<(std::ostream &os, const smat  &m);
template ITPP_EXPORT  std::ostream &operator<<(std::ostream &os, const bmat  &m);

template ITPP_EXPORT  std::istream &operator>>(std::istream &is, mat  &m);
template ITPP_EXPORT  std::istream &operator>>(std::istream &is, cmat &m);
template ITPP_EXPORT std::istream &operator>>(std::istream &is, imat  &m);
template ITPP_EXPORT std::istream &operator>>(std::istream &is, smat  &m);
template ITPP_EXPORT std::istream &operator>>(std::istream &is, bmat  &m);

} // namespace itpp

//! \endcond
