/*!
 * \file 
 * \brief FFT/IFFT test program
 * \author Tony Ottosson, Thomas Eriksson, Simon Wood and Adam Piatyszek
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

#include <itpp/itbase.h>

using namespace itpp;
using namespace std;


#if defined(HAVE_FFTW) || defined(HAVE_MKL)

int main()
{
  int N=pow2i(16), NN=1000, i;
  vec a(N), b(N);
  double *p1=a._data(),*p2=b._data();

  //Remove comments below to measure execution times.
  //Real_Timer time;

  //time.tic();
  for (i=0; i<NN; i++) {
    memcpy(p1,p2,4*N);
  }
  //time.toc();

  N=17, NN=100, i;

  vec in(N), in_final;
  cvec out, out_real, in_complex, inv_out;

  in = randn(N);

  cout << "fft:" << endl;
  //time.tic();
  for (i=0; i<NN; i++) {
    in_complex = to_cvec(in);
    fft(in_complex, out);
  }
  //time.toc();

  cout << "Real fft:" << endl;
  //time.tic();
  for (i=0; i<NN; i++) {
    fft_real(in, out_real);
  }
  //time.toc();

  cout << "Real ifft:" << endl;
  //time.tic();
  for (i=0; i<NN; i++) {
    ifft_real(out_real, in_final);
  }
  //time.toc();

  cout << "Old ifft:" << endl;
  //time.tic();
  for (i=0; i<NN; i++) {
    inv_out = ifft(out);
  }
  //time.toc();

  return 0;
}

#else

int main() { 
  cerr << "Error: FFTW (or MKL) is needed for this test program" << endl; 
  return 1;
}

#endif
