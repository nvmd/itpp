/*!
 * \file
 * \brief Bessel help functions header. For internal use only.
 * \author Tony Ottosson
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

#ifndef BESSEL_INTERNAL_H
#define BESSEL_INTERNAL_H

//! \cond

extern "C" {

  double chbevl(double x, double array[], int n);
  double hyperg(double a, double b, double x);
  int airy(double x, double *ai, double *aip, double *bi, double *bip);
  double polevl(double x, double coef[], int N);
  double p1evl(double x, double coef[], int N);

  double i0(double x);
  double i0e(double x);
  double i1(double x);
  double i1e(double x);

  double k0(double x);
  double k0e(double x);
  double k1(double x);
  double k1e(double x);

  double iv(double nu, double x);
  double jv(double nu, double x);
  double yv(double nu, double x);
  double kn(int n, double x);

} // extern C

//! \endcond

#endif // #ifndef BESSEL_INTERNAL_H
