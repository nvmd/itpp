/*!
 * \file
 * \brief Definitions of Fourier, Hadamard, Walsh-Hadamard, and 2D Hadamard 
 * transforms
 * \author Tony Ottosson, Thomas Eriksson, and Simon Wood
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
 * 
 * This file uses the FFTW package version 3.0.x that is distributed under 
 * the GNU GPL licence. For details see http://www.fftw.org/.
 */

#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>

namespace itpp {

  /*! 
    \addtogroup fft
    \brief One dimensional fast fourier transform
    \author Tony Ottosson

    The functions \code X = fft(x) \endcode and \code x = ifft(X) \endcode are
    the fourier and inverse fourier transforms of size \a N defined as:
    \f[
    X(k) = \sum_{j=0}^{N-1} x(j) e^{-2\pi j k \cdot i / N}
    \f]

    \f[
    x(j) = \frac{1}{N} \sum_{k=0}^{N-1} X(k) e^{2\pi j k \cdot i / N}
    \f]

    \code Y = fft(X, N) \endcode performs zero-padding up to size N and then
    performs an N-size fft.

    The implementation is built upon the FFTW library (www.fftw.org) which is
    a big set of fast and reliable fft routines used e.g. in matlab.

    The routine is fastest for powers of two. Furthermore, the second time you
    call the routine with the same size, the calculation is much faster due to
    many things were calculated and stored the first time the routine was called.
  */

  //!\addtogroup fft
  //!@{

  //! Fast Fourier Transform
  void fft(const cvec &in, cvec &out);
  //! Fast Fourier Transform
  cvec fft(const cvec &in);
  //! Fast Fourier Transform, with zero-padding up to size N
  cvec fft(const cvec &in, const int N);
  //! Inverse Fast Fourier Transform
  void ifft(const cvec &in, cvec &out);
  //! Inverse Fast Fourier Transform
  cvec ifft(const cvec &in);
  //! Inverse Fast Fourier Transform, with zero-padding up to size N
  cvec ifft(const cvec &in, const int N);

  //! Real Fast Fourier Transform
  void fft_real(const vec& in, cvec &out);
  //! Real Fast Fourier Transform
  cvec fft_real(const vec& in);
  //! Real Fast Fourier Transform, with zero-padding up to size N
  cvec fft_real(const vec &in, const int N);
  //! Inverse Real Fast Fourier Transform. Assumes even size.
  void ifft_real(const cvec &in, vec &out);
  //! Inverse Real Fast Fourier Transform. Assumes even size.
  vec ifft_real(const cvec &in);
  //! Inverse Real Fast Fourier Transform, with zero-padding up to size N
  vec ifft_real(const cvec &in, const int N);
  //!@}


  /*! 
		\addtogroup dct
    \brief One dimensional Dicrete Cosine Transform
    \author Tony Ottosson

    The functions \code X = dct(x) \endcode and \code x = idct(X) \endcode are
    the dicrete cosine and inverse discrete cosine transforms of size \a N defined as:
    \f[
    X(k) = w(k) \sum_{j=0}^{N-1} x(j) \cos \left( \frac{(2j+1) k \pi}{2 N} \right)
    \f]

    \f[
    x(j) = \sum_{k=0}^{N-1} w(k) X(k) \cos \left( \frac{(2j+1) k \pi}{2 N} \right)
    \f]
    where \f$w(k) = 1/sqrt{N}\f$ for \f$k=0\f$ and \f$w(k) = sqrt{2/N}\f$ for \f$k\geq 1\f$.


    The implementation is built upon the FFTW library (www.fftw.org) which is
    a big set of fast and reliable fft routines used e.g. in matlab.

    The routine is fastest for powers of two. Furthermore, the second time you
    call the routine with the same size, the calculation is much faster due to
    many things were calculated and stored the first time the routine was called.
  */

  //!\addtogroup dct
  //!@{
  //! Discrete Cosine Transform (DCT)
  void dct(const vec &in, vec &out);
  //! Discrete Cosine Transform (DCT)
  vec dct(const vec &in);
  //! Inverse Discrete Cosine Transform (IDCT)
  void idct(const vec &in, vec &out);
  //! Inverse Discrete Cosine Transform (IDCT)
  vec idct(const vec &in);
  //!@}
  
	
  //!\addtogroup fht
  //!@{

  //! Fast Hadamard Transform
  template <class T> Vec<T> dht(const Vec<T> &v);
  //! Fast Hadamard Transform
  template <class T> void dht(const Vec<T> &vin, Vec<T> &vout);
  //! Fast Hadamard Transform - memory efficient. Stores the result in \c v
  template <class T> void self_dht(Vec<T> &v);

  //! Fast Walsh Hadamard Transform
  template <class T> Vec<T> dwht(const Vec<T> &v);
  //! Fast Walsh Hadamard Transform
  template <class T> void dwht(const Vec<T> &vin, Vec<T> &vout);
  //! Fast Walsh Hadamard Transform - memory efficient. Stores the result in \c v
  template <class T> void self_dwht(Vec<T> &v);

  //! Fast 2D Hadamard Transform
  template <class T> Mat<T> dht2(const Mat<T> &m);
  //! Fast 2D Walsh Hadamard Transform
  template <class T> Mat<T> dwht2(const Mat<T> &m);
  //!@}

  /////////////////////////////
  // template instantiation
  /////////////////////////////
#ifndef _MSC_VER

  //! Template instantiation of dht
  extern template vec dht(const vec &v);
  //! Template instantiation of dht
  extern template cvec dht(const cvec &v);

  //! Template instantiation of dht
  extern template void dht(const vec &vin, vec &vout);
  //! Template instantiation of dht
  extern template void dht(const cvec &vin, cvec &vout);

  //! Template instantiation of self_dht
  extern template void self_dht(vec &v);
  //! Template instantiation of self_dht
  extern template void self_dht(cvec &v);

  //! Template instantiation of dwht
  extern template vec dwht(const vec &v);
  //! Template instantiation of dwht
  extern template cvec dwht(const cvec &v);

  //! Template instantiation of dwht
  extern template void dwht(const vec &vin, vec &vout);
  //! Template instantiation of dwht
  extern template void dwht(const cvec &vin, cvec &vout);

  //! Template instantiation of self_dwht
  extern template void self_dwht(vec &v);
  //! Template instantiation of self_dwht
  extern template void self_dwht(cvec &v);

  //! Template instantiation of dht2
  extern template mat  dht2(const mat &m);
  //! Template instantiation of dht2
  extern template cmat dht2(const cmat &m);

  //! Template instantiation of dht2
  extern template mat  dwht2(const mat &m);
  //! Template instantiation of dht2
  extern template cmat dwht2(const cmat &m);

#endif // #ifndef _MSC_VER

} // namespace itpp

#endif // #ifndef TRANSFORMS_H
