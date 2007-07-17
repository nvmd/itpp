/*!
 * \file
 * \brief Implementation of window functions
 * \author Tony Ottosson, Tobias Ringstrom, Pal Frenger and Adam Piatyszek
 *
 * $Date$
 * $Revision$
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

#include <itpp/signal/window.h>


namespace itpp {


  vec hamming(int n)
  {
    vec	t(n);

    if (n == 1)
      t(0) = 0.08;
    else
      for (int i=0;i<n;i++)
	t[i]=(0.54-0.46*std::cos(2.0*pi*i/(n-1)));

    return t;
  }

  vec hanning(int n)
  {
    vec	t(n);

    for (int i=0;i<n;i++)
      t(i) = 0.5 * (1.0 - std::cos(2.0*pi*(i+1)/(n+1)));

    return t;
  }

  // matlab version
  vec hann(int n)
  {
    vec	t(n);

    for (int i=0;i<n;i++)
      t(i) = 0.5 * (1.0 - std::cos(2.0*pi*i/(n-1)));

    return t;
  }

  vec blackman(int n)
  {
    vec	t(n);

    for (int i=0;i<n;i++)
      t(i) = 0.42 - 0.5 * std::cos(2.0*pi*i/(n-1)) + 0.08 * std::cos(4.0*pi*i/(n-1));

    return t;
  }

  vec triang(int n)
  {
    vec	t(n);

    if (n % 2) { // Odd
      for (int i=0; i<n/2; i++)
	t(i) = t(n-i-1) = 2.0*(i+1)/(n+1);
      t(n/2) = 1.0;
    } else
      for (int i=0; i<n/2; i++)
	t(i) = t(n-i-1) = (2.0*i+1)/n;

    return t;
  }

  vec sqrt_win(int n)
  {
    vec	t(n);

    if (n % 2) { // Odd
      for (int i=0; i<n/2; i++)
	t(i) = t(n-i-1) = std::sqrt(2.0*(i+1)/(n+1));
      t(n/2) = 1.0;
    } else
      for (int i=0; i<n/2; i++)
	t(i) = t(n-i-1) = std::sqrt((2.0*i+1)/n);

    return t;
  }



} // namespace itpp
