/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2003 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Bessel help functions header. For internal use only.
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __bessel_internal_h
#define __bessel_internal_h


#ifndef DOXYGEN_SHOULD_SKIP_THIS

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


#endif //DOXYGEN_SHOULD_SKIP_THIS

#endif // __bessel_internal_h

