/*!
 * \file
 * \brief Definitions of window functions
 * \author Tony Ottosson, Tobias Ringstrom, Pal Frenger and Adam Piatyszek
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

#ifndef WINDOW_H
#define WINDOW_H

#include <itpp/base/vec.h>


namespace itpp {

  /*!
    \addtogroup windfunc
  */

  /*!\addtogroup windfunc
    \brief Windowing functions

  */
  //!@{

  /*! \brief Hamming window

  The \c n size Hamming window is a vector \f$w\f$ where the \f$i\f$th component is
  \f[
  w_i = 0.54 - 0.46 \cos(2\pi i/(n-1))
  \f]
  */
  vec hamming(int size);


  /*! \brief Hanning window

  The \c n size Hanning window is a vector \f$w\f$ where the \f$i\f$th component is
  \f[
  w_i = 0.5(1 - \cos(2\pi (i+1)/(n+1))
  \f]

  Observe that this function is not the same as the hann() function which is defined
  as in matlab.
  */
  vec hanning(int n);

  /*! \brief Hanning window compatible with matlab

  The \c n size Hanning window is a vector \f$w\f$ where the \f$i\f$th component is
  \f[
  w_i = 0.5(1 - \cos(2\pi i/(n-1))
  \f]
  */
  vec hann(int n);


  /*! \brief Blackman window

  The \c n size Blackman window is a vector \f$w\f$ where the \f$i\f$th component is
  \f[
  w_i = 0.42 - 0.5\cos(2\pi i/(n-1)) + 0.08\cos(4\pi i/(n-1))
  \f]
  */
  vec blackman(int n);

  /*! \brief Triangular window


  The \c n size triangle window is a vector \f$w\f$ where the \f$i\f$th component is
  \f[
  w_i = w_{n-i-1} = \frac{2(i+1)}{n+1}
  \f]
  for \c n odd and for \c n even
  \f[
  w_i = w_{n-i-1} = \frac{2i+1}{n}
  \f]
  */
  vec triang(int n);

  /*! \brief Square root window

  The square-root of the Triangle window.
  sqrt_win(n) = sqrt(triang(n))
  */
  vec sqrt_win(int n);
  //!@}


} //namespace itpp

#endif // #ifndef WINDOW_H
