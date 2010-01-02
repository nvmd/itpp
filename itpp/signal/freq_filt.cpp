/*!
 * \file
 * \brief Implementation of Freq_Filt Class
 * \author Simon Wood
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

#include <itpp/signal/freq_filt.h>
#include <itpp/signal/transforms.h>
#include <itpp/base/math/elem_math.h>

//! \cond

namespace itpp
{

// Overlap-add routine
template<class Num_T>
void Freq_Filt<Num_T>::overlap_add(const cvec&x, cvec &y)
{
  int nb = impulse.length();
  int nx = x.length();

  y.set_size(nx, false);
  y.zeros();
  cvec X, Y;
  int istart = 0;
  int L = blksize;
  while (istart < nx) {
    int iend = std::min(istart + L - 1, nx - 1);

    X = fft(x(istart, iend), fftsize);
    Y = ifft(elem_mult(X, B));
    Y.set_subvector(0, Y(0, nb - 2) + zfinal);
    int yend = std::min(nx - 1, istart + fftsize - 1);
    y.set_subvector(istart, Y(0, yend - istart));
    zfinal = Y(fftsize - (nb - 1), fftsize - 1);
    istart += L;
  }
}

template<>
vec Freq_Filt<double>::overlap_add(const vec &x)
{
  cvec y; // Size of y is set later
  overlap_add(to_cvec(x), y);
  return real(y);
}

template<>
svec Freq_Filt<short>::overlap_add(const svec &x)
{
  cvec y; // Size of y is set later
  overlap_add(to_cvec(x), y);
  return to_svec(real(y));
}

template<>
ivec Freq_Filt<int>::overlap_add(const ivec &x)
{
  cvec y; // Size of y is set later
  overlap_add(to_cvec(x), y);
  return to_ivec(real(y));
}

template<>
cvec Freq_Filt<std::complex<double> >::overlap_add(const cvec &x)
{
  cvec y; // Size of y is set later
  overlap_add(x, y);
  return y;
}

} // namespace itpp

//! \endcond
