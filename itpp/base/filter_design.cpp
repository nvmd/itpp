/*!
 * \file
 * \brief Filter design functions
 * \author Tony Ottosson
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

#include <itpp/base/filter_design.h>

#include <itpp/base/poly.h>
#include <itpp/base/elmatfunc.h>
#include <itpp/base/transforms.h>


namespace itpp {


  void polystab(const vec &a, vec &out)
  {
    cvec r;
    roots(a, r);
    
    for (int i=0; i<r.size(); i++) {
      if (abs(r(i)) > 1)
	r(i) = std::complex<double>(1.0)/conj(r(i));
    }
    out = real(std::complex<double>(a(0)) * poly(r));
  }

  void polystab(const cvec &a, cvec &out)
  {
    cvec r;
    roots(a, r);
    
    for (int i=0; i<r.size(); i++) {
      if (abs(r(i)) > 1)
	r(i) = std::complex<double>(1.0)/conj(r(i));
    }
    out = a(0) * poly(r);
  }

  void freqz(const cvec &b, const cvec &a, const int N, cvec &h, vec &w)
  {
    w = pi*linspace(0, N-1, N)/double(N);

    cvec ha, hb;
    hb = fft( b, 2*N );
    hb = hb(0, N-1);

    ha = fft( a, 2*N );
    ha = ha(0, N-1);

    h = elem_div(hb, ha);
  }

  cvec freqz(const cvec &b, const cvec &a, const int N)
  {
    cvec h;
    vec w;

    freqz(b, a, N, h, w);

    return h;
  }


  cvec freqz(const cvec &b, const cvec &a, const vec &w)
  {
    int la = a.size(), lb = b.size(), k = std::max(la, lb);

    cvec h, ha, hb;

    // Evaluate the nominator and denominator at the given frequencies
    hb = polyval( zero_pad(b, k), to_cvec(cos(w), sin(w)) );
    ha = polyval( zero_pad(a, k), to_cvec(cos(w), sin(w)) );

    h = elem_div(hb, ha);

    return h;
  }


} // namespace itpp
