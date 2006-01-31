/*!
* \file 
* \brief Matrix inversion routines test program
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

#include <itpp/itbase.h>

using namespace itpp;
using namespace std;


#if defined(HAVE_LAPACK)

int main(void)
{
 cout << "======================================" << endl;
 cout << "    Test of Matrix inversion routines" << endl;
 cout << "======================================" << endl;

 {
   cout << "Real matrix" << endl;

   mat X = randn(5, 5), Y;
   Y = inv(X);
   cout << "X = " << round_to_zero(X) << endl;
   cout << "inv(X) = " << round_to_zero(Y) << endl;

   X = randn(5, 5);
   Y = inv(X);
   cout << "X = " << round_to_zero(X) << endl;
   cout << "inv(X) = " << round_to_zero(Y) << endl;
 }
 {
   cout << endl << "Complex matrix" << endl;

   cmat X = randn_c(5, 5), Y;
   Y = inv(X);
   cout << "X = " << round_to_zero(X) << endl;
   cout << "inv(X) = " << round_to_zero(Y) << endl;

   X = randn_c(5, 5);
   Y = inv(X);
   cout << "X = " << round_to_zero(X) << endl;
   cout << "inv(X) = " << round_to_zero(Y) << endl;
 }

 return 0;
}

#else

int main() { 
 cerr << "Error: LAPACK library is needed to run this test program" << endl; 
 return 1;
}

#endif
