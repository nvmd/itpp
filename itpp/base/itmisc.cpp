/*!
 * \file
 * \brief Implementation of IT++ miscleaneous functions
 * \author Tony Ottosson and Adam Piatyszek
 * 
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
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
 *
 * -------------------------------------------------------------------------
 */

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <itpp/base/itmisc.h>
#include <itpp/base/vec.h>


namespace itpp { 

  //Calculates factorial coefficient for index <= 170.
  double fact(int index)
  {
    it_error_if(index > 170, "\nThe function double factfp(int index) overflows if index > 170. \nUse your head instead!");
    it_error_if(index < 0, "\nThe function double factfp(int index) cannot evaluate if index. < 0");
    double prod = 1;
    for (int i=1; i<=index; i++)
      prod *= double(i);
    return prod;
  }


  long gcd(long a, long b)
  {
    long v, u, t, q;

    it_assert(a>=0,"long gcd(long a, long b): a and b must be non-negative integers");
    it_assert(b>=0,"long gcd(long a, long b): a and b must be non-negative integers");

    u = std::abs(a);
    v = std::abs(b);
    while (v>0) {
      q = u / v;
      t = u - v*q;
      u = v;
      v = t;
    }
    return(u);
  }


  // Calculates binomial coefficient "n over k".
  double binom(int n, int k) {
    it_error_if(k>n,"Error in double binom(int n, int k).\nn must be larger than k.");
    k = n-k<k ? n-k : k;

    vec talj(k), namn(k);
    int i;
    double out = 1.0;
    for (i=0; i<k; i++)
      namn(i) = double(i+1);

    int pos = 0;
    for (i=n; i>=(n-k+1); i--) {
      talj(pos) = double(i);
      pos++;
    }

    for (i=0; i<k; i++) {
      out *= talj(i) / namn(k-1-i);
    }
    return ( out );
  }

  int binom_i(int n, int k)
  {
    ivec v(n);
    int i, j;

    if (n > (k+1)/2)
      n = k+1-n;

    v = 0;
    v(0) = 1;

    for (i=0; i<k-n; i++)
      for (j=1; j<n; j++)
	v(j) += v(j-1);

    return v(n-1);
  }

  // Calculates the base 10-logarithm of the binomial coefficient "n over k".
  double log_binom(int n, int k) {
    it_error_if(k>n,"Error in double log_binom(int n, int k).\nn must be larger than k.");
    k = n-k<k ? n-k : k;

    int i;
    double out = 0.0;
    for (i=0; i<k; i++)
      out += log10((double)(n-i)) - log10((double)(i+1));

    return out;
  }


  std::string itpp_version(void)
  {
#ifdef PACKAGE_VERSION
    return std::string(PACKAGE_VERSION);
#else
    return std::string("Warning: Version unknown!");
#endif
  }

  bool check_big_endianness()
  {
    int i = 1;
    char *p = (char *) &i;
    if (p[0] == 1) // Lowest address contains the least significant byte
      return false; // LITTLE_ENDIAN
    else
      return true; // BIG_ENDIAN
  }
  
} //namespace itpp
