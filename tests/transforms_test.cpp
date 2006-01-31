/*!
* \file 
* \brief Transforms test program
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


#if defined(HAVE_FFT)

int main()
{
 cout << "========================" << endl;
 cout << "   Test of Transforms   " << endl;
 cout << "========================" << endl;

 {
   cout << "Real vector: fft_real(x,y), ifft_real(y,z)" << endl;
   int N = 16;

   vec x, z;
   cvec y;
   
   x = randn(N);
   fft_real(x, y);
   ifft_real(y, z);
   
   cout << "x = " << x << endl;
   cout << "y = " << y << endl;
   cout << "z = " << z << endl;
 }
 {
   cout << "Complex vector: fft(x,y), ifft(y,z)" << endl;
   int N = 16;

   cvec x, y, z;
   
   x = randn_c(N);
   fft(x, y);
   ifft(y, z);
   
   cout << "x = " << x << endl;
   cout << "y = " << y << endl;
   cout << "z = " << z << endl;
 }

 return 0;
}

#else

int main() { 
 cerr << "Error: FFTW library is needed to run this test program" << endl; 
 return 1;
}

#endif
