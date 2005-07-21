/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2005 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Freq_Filt Class Implementation
  \author Simon Wood

  $Revision$

  $Date$ 
*/

#include <itpp/base/freq_filt.h>
#include <itpp/base/transforms.h>
#include <itpp/base/elmatfunc.h>

namespace itpp {

  // Overlap-add routine
  template<class Num_T>
  void Freq_Filt<Num_T>::overlap_add(const cvec&x, cvec &y)
  {
    int nb = impulse.length();
    int nx = x.length();

    y.set_size(nx,false);
    y.zeros();
    cvec X,Y;
    int istart = 0;
    int L = blksize;
    while(istart < nx)
      {
	int iend = std::min(istart+L-1,nx-1);

	X = fft(x(istart,iend),fftsize);
	Y = ifft(elem_mult(X,B));
	Y.set_subvector(0,nb-2,Y(0,nb-2) + zfinal);
	int yend = std::min(nx-1,istart+fftsize-1);
	y.set_subvector(istart,yend,Y(0,yend-istart));
	zfinal = Y(fftsize-(nb-1),fftsize-1);
	istart += L;
      }
  }

  template<>
  vec Freq_Filt<double>::overlap_add(const vec &x)
  {
    cvec y; // Size of y is set later
    overlap_add(to_cvec(x),y);
    return real(y);
  }

  template<>
  svec Freq_Filt<short>::overlap_add(const svec &x)
  {
    cvec y; // Size of y is set later
    overlap_add(to_cvec(x),y);
    return to_svec(real(y));
  }

  template<>
  ivec Freq_Filt<int>::overlap_add(const ivec &x)
  {
    cvec y; // Size of y is set later
    overlap_add(to_cvec(x),y);
    return to_ivec(real(y));
  }

  template<>
  cvec Freq_Filt<std::complex<double> >::overlap_add(const cvec &x)
  {
    cvec y; // Size of y is set later
    overlap_add(x,y);
    return y;
  }

} // namespace itpp
