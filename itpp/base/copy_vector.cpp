/*!
 * \file
 * \brief Vector copy functions for internal use
 * \author Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2008  (see AUTHORS file for a list of contributors)
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

#include <itpp/base/copy_vector.h>

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#if defined (HAVE_BLAS)
#  include <itpp/base/blas.h>
#endif


namespace itpp
{

void copy_vector(int n, const double *x, double *y)
{
#if defined (HAVE_BLAS)
  int_blas_t incr = 1;
  int_blas_t _n = static_cast<int_blas_t>(n);
  blas::dcopy_(&_n, x, &incr, y, &incr);
#else
  memcpy(y, x, static_cast<unsigned int>(n) * sizeof(double));
#endif
}

void copy_vector(int n, const std::complex<double> *x, std::complex<double> *y)
{
#if defined (HAVE_BLAS)
  int_blas_t incr = 1;
  int_blas_t _n = static_cast<int_blas_t>(n);
  blas::zcopy_(&_n, x, &incr, y, &incr);
#else
  memcpy(y, x, static_cast<unsigned int>(n) * sizeof(std::complex<double>));
#endif
}

void copy_vector(int n, const double *x, int incx, double *y, int incy)
{
#if defined (HAVE_BLAS)
  int_blas_t _n = static_cast<int_blas_t>(n);
  int_blas_t _incx = static_cast<int_blas_t>(incx);
  int_blas_t _incy = static_cast<int_blas_t>(incy);
  blas::dcopy_(&_n, x, &_incx, y, &_incy);
#else
  for (int i = 0; i < n; i++)
    y[i*incy] = x[i*incx];
#endif
}

void copy_vector(int n, const std::complex<double> *x, int incx,
                 std::complex<double> *y, int incy)
{
#if defined (HAVE_BLAS)
  int_blas_t _n = static_cast<int_blas_t>(n);
  int_blas_t _incx = static_cast<int_blas_t>(incx);
  int_blas_t _incy = static_cast<int_blas_t>(incy);
  blas::zcopy_(&_n, x, &_incx, y, &_incy);
#else
  for (int i = 0; i < n; i++)
    y[i*incy] = x[i*incx];
#endif
}

void swap_vector(int n, double *x, double *y)
{
#if defined (HAVE_BLAS)
  int_blas_t incr = 1;
  int_blas_t _n = static_cast<int_blas_t>(n);
  blas::dswap_(&_n, x, &incr, y, &incr);
#else
  for (int i = 0; i < n; i++)
    std::swap(x[i], y[i]);
#endif
}

void swap_vector(int n, std::complex<double> *x, std::complex<double> *y)
{
#if defined (HAVE_BLAS)
  int_blas_t incr = 1;
  int_blas_t _n = static_cast<int_blas_t>(n);
  blas::zswap_(&_n, x, &incr, y, &incr);
#else
  for (int i = 0; i < n; i++)
    std::swap(x[i], y[i]);
#endif
}

void swap_vector(int n, double *x, int incx, double *y, int incy)
{
#if defined (HAVE_BLAS)
  int_blas_t _n = static_cast<int_blas_t>(n);
  int_blas_t _incx = static_cast<int_blas_t>(incx);
  int_blas_t _incy = static_cast<int_blas_t>(incy);
  blas::dswap_(&_n, x, &_incx, y, &_incy);
#else
  for (int i = 0; i < n; i++)
    std::swap(x[i*incx], y[i*incy]);
#endif
}

void swap_vector(int n, std::complex<double> *x, int incx,
                 std::complex<double> *y, int incy)
{
#if defined (HAVE_BLAS)
  int_blas_t _n = static_cast<int_blas_t>(n);
  int_blas_t _incx = static_cast<int_blas_t>(incx);
  int_blas_t _incy = static_cast<int_blas_t>(incy);
  blas::zswap_(&_n, x, &_incx, y, &_incy);
#else
  for (int i = 0; i < n; i++)
    std::swap(x[i*incx], y[i*incy]);
#endif
}

void scal_vector(int n, double alpha, double *x)
{
#if defined(HAVE_BLAS)
  int_blas_t incr = 1;
  int_blas_t _n = static_cast<int_blas_t>(n);
  blas::dscal_(&_n, &alpha, x, &incr);
#else
  if (alpha != 1.0) {
    for (int i = 0; i < n; ++i) {
      x[i] *= alpha;
    }
  }
#endif
}

void scal_vector(int n, std::complex<double> alpha, std::complex<double> *x)
{
#if defined(HAVE_BLAS)
  int_blas_t incr = 1;
  int_blas_t _n = static_cast<int_blas_t>(n);
  blas::zscal_(&_n, &alpha, x, &incr);
#else
  if (alpha != std::complex<double>(1.0)) {
    for (int i = 0; i < n; ++i) {
      x[i] *= alpha;
    }
  }
#endif
}

void scal_vector(int n, double alpha, double *x, int incx)
{
#if defined(HAVE_BLAS)
  int_blas_t _n = static_cast<int_blas_t>(n);
  int_blas_t _incx = static_cast<int_blas_t>(incx);
  blas::dscal_(&_n, &alpha, x, &_incx);
#else
  if (alpha != 1.0) {
    for (int i = 0; i < n; ++i) {
      x[i*incx] *= alpha;
    }
  }
#endif
}

void scal_vector(int n, std::complex<double> alpha, std::complex<double> *x,
                 int incx)
{
#if defined(HAVE_BLAS)
  int_blas_t _n = static_cast<int_blas_t>(n);
  int_blas_t _incx = static_cast<int_blas_t>(incx);
  blas::zscal_(&_n, &alpha, x, &_incx);
#else
  if (alpha != std::complex<double>(1.0)) {
    for (int i = 0; i < n; ++i) {
      x[i*incx] *= alpha;
    }
  }
#endif
}

void axpy_vector(int n, double alpha, const double *x, double *y)
{
#if defined(HAVE_BLAS)
  int_blas_t incr = 1;
  int_blas_t _n = static_cast<int_blas_t>(n);
  blas::daxpy_(&_n, &alpha, x, &incr, y, &incr);
#else
  if (alpha != 1.0) {
    for (int i = 0; i < n; ++i) {
      y[i] += alpha * x[i];
    }
  }
  else {
    for (int i = 0; i < n; ++i) {
      y[i] += x[i];
    }
  }
#endif
}

void axpy_vector(int n, std::complex<double> alpha,
                 const std::complex<double> *x, std::complex<double> *y)
{
#if defined(HAVE_BLAS)
  int_blas_t incr = 1;
  int_blas_t _n = static_cast<int_blas_t>(n);
  blas::zaxpy_(&_n, &alpha, x, &incr, y, &incr);
#else
  if (alpha != std::complex<double>(1.0)) {
    for (int i = 0; i < n; ++i) {
      y[i] += alpha * x[i];
    }
  }
  else {
    for (int i = 0; i < n; ++i) {
      y[i] += x[i];
    }
  }
#endif
}

void axpy_vector(int n, double alpha, const double *x, int incx, double *y,
                 int incy)
{
#if defined(HAVE_BLAS)
  int_blas_t _n = static_cast<int_blas_t>(n);
  int_blas_t _incx = static_cast<int_blas_t>(incx);
  int_blas_t _incy = static_cast<int_blas_t>(incy);
  blas::daxpy_(&_n, &alpha, x, &_incx, y, &_incy);
#else
  if (alpha != 1.0) {
    for (int i = 0; i < n; ++i) {
      y[i*incy] += alpha * x[i*incx];
    }
  }
  else {
    for (int i = 0; i < n; ++i) {
      y[i*incy] += x[i*incx];
    }
  }
#endif
}

void axpy_vector(int n, std::complex<double> alpha,
                 const std::complex<double> *x, int incx,
                 std::complex<double> *y, int incy)
{
#if defined(HAVE_BLAS)
  int_blas_t _n = static_cast<int_blas_t>(n);
  int_blas_t _incx = static_cast<int_blas_t>(incx);
  int_blas_t _incy = static_cast<int_blas_t>(incy);
  blas::zaxpy_(&_n, &alpha, x, &_incx, y, &_incy);
#else
  if (alpha != std::complex<double>(1.0)) {
    for (int i = 0; i < n; ++i) {
      y[i*incy] += alpha * x[i*incx];
    }
  }
  else {
    for (int i = 0; i < n; ++i) {
      y[i*incy] += x[i*incx];
    }
  }
#endif
}

} // namespace itpp
