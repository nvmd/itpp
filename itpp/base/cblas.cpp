/*!
 * \file
 * \brief C interface for a limited set of BLAS functions from ACML
 * \author Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

#include <itpp/base/cblas.h>

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#ifdef HAVE_CBLAS_ACML
#include <acml.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

  void cblas_dcopy(const int N, const double *X, const int incX,
		   double *Y, const int incY)
  {
    dcopy(N, (double*)(X), incX, Y, incY);
  }

  void cblas_zcopy(const int N, const void *X, const int incX,
		   void *Y, const int incY)
  {
    zcopy(N, (doublecomplex*)(X), incX, (doublecomplex*)(Y), incY);
  }

  void cblas_dswap(const int N, double *X, const int incX,
		   double *Y, const int incY)
  {
    dswap(N, X, incX, Y, incY);
  }

  void cblas_zswap(const int N, void *X, const int incX,
		   void *Y, const int incY)
  {
    zswap(N, (doublecomplex*)(X), incX, (doublecomplex*)(Y), incY);
  }

  double cblas_ddot(const int N, const double *X, const int incX,
		    const double *Y, const int incY)
  {
    return ddot(N, (double*)(X), incX, (double*)(Y), incY);
  }

  void cblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
  {
    doublecomplex r;
    r = zdotu(N, (doublecomplex*)(X), incX, (doublecomplex*)(Y), incY);
    *((double*)(dotu)) = r.real;
    *((double*)(dotu) + 1) = r.imag;
  }

  void cblas_dgemm(const enum CBLAS_ORDER Order,
		   const enum CBLAS_TRANSPOSE TransA,
		   const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
		   const int K, const double alpha, const double *A,
		   const int lda, const double *B, const int ldb,
		   const double beta, double *C, const int ldc)
  {
    char trans_a, trans_b;
    switch (TransA) {
    case CblasNoTrans: trans_a = 'n'; break;
    case CblasTrans: trans_a = 't'; break;
    case CblasConjTrans: trans_a = 'c'; break;
    case AtlasConj: trans_a = 'c'; break;
    }
    switch (TransB) {
    case CblasNoTrans: trans_b = 'n'; break;
    case CblasTrans: trans_b = 't'; break;
    case CblasConjTrans: trans_b = 'c'; break;
    case AtlasConj: trans_b = 'c'; break;
    }
    dgemm(trans_a, trans_b, M, N, K, alpha, (double*)(A), lda, (double*)(B),
	  ldb, beta, C, ldc);
  }

  void cblas_zgemm(const enum CBLAS_ORDER Order,
		   const enum CBLAS_TRANSPOSE TransA,
		   const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
		   const int K, const void *alpha, const void *A,
		   const int lda, const void *B, const int ldb,
		   const void *beta, void *C, const int ldc)
  {
    char trans_a, trans_b;
    switch (TransA) {
    case CblasNoTrans: trans_a = 'n'; break;
    case CblasTrans: trans_a = 't'; break;
    case CblasConjTrans: trans_a = 'c'; break;
    case AtlasConj: trans_a = 'c'; break;
    }
    switch (TransB) {
    case CblasNoTrans: trans_b = 'n'; break;
    case CblasTrans: trans_b = 't'; break;
    case CblasConjTrans: trans_b = 'c'; break;
    case AtlasConj: trans_b = 'c'; break;
    }
    zgemm(trans_a, trans_b, M, N, K, (doublecomplex*)(alpha),
	  (doublecomplex*)(A), lda, (doublecomplex*)(B), ldb,
	  (doublecomplex*)(beta), (doublecomplex*)(C), ldc);
  }

  void cblas_dgemv(const enum CBLAS_ORDER Order,
		   const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
		   const double alpha, const double *A, const int lda,
		   const double *X, const int incX, const double beta,
		   double *Y, const int incY)
  {
    char trans_a;
    switch (TransA) {
    case CblasNoTrans: trans_a = 'n'; break;
    case CblasTrans: trans_a = 't'; break;
    case CblasConjTrans: trans_a = 'c'; break;
    case AtlasConj: trans_a = 'c'; break;
    }
    dgemv(trans_a, M, N, alpha, (double*)(A), lda, (double*)(X), incX, beta,
	  Y, incY);
  }

  void cblas_zgemv(const enum CBLAS_ORDER Order,
		   const enum CBLAS_TRANSPOSE TransA, const int M, const int N,
		   const void *alpha, const void *A, const int lda,
		   const void *X, const int incX, const void *beta,
		   void *Y, const int incY)
  {
    char trans_a;
    switch (TransA) {
    case CblasNoTrans: trans_a = 'n'; break;
    case CblasTrans: trans_a = 't'; break;
    case CblasConjTrans: trans_a = 'c'; break;
    case AtlasConj: trans_a = 'c'; break;
    }
    zgemv(trans_a, M, N, (doublecomplex*)(alpha), (doublecomplex*)(A), lda,
	  (doublecomplex*)(X), incX, (doublecomplex*)(beta),
	  (doublecomplex*)(Y), incY);
  }

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* #ifdef HAVE_CBLAS_ACML */
