/*!
* \file 
* \brief Statistical routines test program
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

int main()
{
 cout << "=================================" << endl;
 cout << "  Test of statistical routines" << endl;
 cout << "=================================" << endl;

 vec a = randn(5);

 cout << "a = " << a << endl;
 cout << "max(a) = " << max(a) << endl;
 cout << "min(a) = " << min(a) << endl;

 mat A = randn(5,5);
 cout << "A = " << A << endl << endl;

 cout << "max(A) = " << max(A) << endl;
 cout << "max(A,1) = " << max(A,1) << endl;
 cout << "max(A,2) = " << max(A,2) << endl;
 cout << "min(A) = " << min(A) << endl;
 cout << "min(A,1) = " << min(A,1) << endl;
 cout << "min(A,2) = " << min(A,2) << endl;

 cout << "norm(A) = " << norm(A) << endl;
 cout << "norm(A,2) = " << norm(A,2) << endl;
 cout << "norm(A,1) = " << norm(A,1) << endl;

 return 0;
}

#else

int main() { 
 cerr << "Error: LAPACK library is needed to run this test program" << endl; 
 return 1;
}

#endif
