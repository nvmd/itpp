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

#ifndef FILTER_DESIGN_H
#define FILTER_DESIGN_H

#include <complex>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include <itpp/itconfig.h>
#include <itpp/base/vec.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/specmat.h>
#include <itpp/base/itassert.h>

namespace itpp {

  /*!
    \addtogroup filters
  */



  /*!
    \brief Polynomial Stabilization
    \ingroup filters
    \author Tony Ottosson

    Stabilizes the polynomial transfer function by replacing all roots outside
    the unit cirlce with their reflection inside the unit circle.
  */
  //!@{
  void polystab(const vec &a, vec &out);
  inline vec polystab(const vec &a) { vec temp; polystab(a, temp); return temp; }
  void polystab(const cvec &a, cvec &out);
  inline cvec polystab(const cvec &a) { cvec temp; polystab(a, temp); return temp; }
  //!@}
  
  /*!
    \brief Frequency response of filter
    \ingroup filters
    \author Tony Ottosson
    
    Calculates the N-point frequency response of the supplied digital filter over the frequencies w.
    If w is not given the response is evaluated over the range 0 to \f$\pi\f$ with N values.
    The default value of N is 512.

    If \c w is supplied polyval() is used. Otherwise the calculation is based on the fft.
  */
  //!@{
  void freqz(const cvec &b, const cvec& a, const int N, cvec &h, vec &w);
  cvec freqz(const cvec &b, const cvec& a, const int N = 512);
  cvec freqz(const cvec &b, const cvec& a, const vec &w);
  //!@}


} // namespace itpp

#endif // #ifndef FILTER_DESIGN_H
