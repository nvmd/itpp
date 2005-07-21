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
  \brief Implementation of Bessel functions.
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#include <cmath>
//#include <iostream> //This line is needed for the Sun WorkShop Compiler 5.0 (CC). Dunno why... (PF 2001-12-03)

#include <itpp/base/bessel.h>
#include <itpp/base/bessel/bessel_internal.h>

namespace itpp { 

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

  // Bessel function of order nu
  double besselj(int nu, double x) { return jn(nu, x); }

  vec besselj(int nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
      out(i) = jn(nu, x(i));

    return out;
  }

  // Bessel function of order nu. nu is real.
  double besselj(double nu, double x) { return jv(nu, x); }

  vec besselj(double nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
      out(i) = jv(nu, x(i));

    return out;
  }

  // Bessel function of second kind of order nu
  double bessely(int nu, double x) { return yn(nu, x); }

  vec bessely(int nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
      out(i) = yn(nu, x(i));

    return out;
  }
  // Bessel function of second kind of order nu
  double bessely(double nu, double x) { return yv(nu, x); }

  vec bessely(double nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
      out(i) = yv(nu, x(i));

    return out;
  }

  // Modified Bessel function of order nu
  double besseli(double nu, double x) { return iv(nu, x); }

  vec besseli(double nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
      out(i) = iv(nu, x(i));

    return out;
  }

  // Modified Bessel function of second kind of order n
  double besselk(int n, double x) { return kn(n, x); }

  vec besselk(int nu, const vec &x)
  {
    vec out(x.size());
    for (int i=0; i<x.size(); i++)
      out(i) = kn(nu, x(i));

    return out;
  }

} //namespace itpp
