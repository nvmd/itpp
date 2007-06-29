/*!
* \file 
* \brief Polynomial routines test program
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

#include <itpp/itbase.h>
#include <itpp/itsignal.h>

using namespace itpp;
using namespace std;


#if defined(HAVE_LAPACK)

int main()
{

 cout << "===================================" << endl;
 cout << "    Test of polynomial routines    " << endl;
 cout << "===================================" << endl;

 {
   cout << "Real polynomials" << endl;
   vec r = randn(3);
   vec p = poly(r);
   cvec r2 = roots(p);
   cout << "Roots, r = " << r << endl;
   cout << "Polynomial, p = " << p << endl;
   cout << "r = roots(p) = " << r2 << endl;

   r = randn(7);
   p = poly(r);
   r2 = roots(p);
   cout << "Roots, r = " << r << endl;
   cout << "Polynomial, p = " << p << endl;
   cout << "r = roots(p) = " << r2 << endl;

   vec x = randn(9);
   vec y = polyval(p, x);
   cout << "x = " << x << endl;
   cout << "polyval(p, x) = " << y << endl;
 }

 {
   cout << "Complex polynomials" << endl;
   cvec r = randn_c(3);
   cvec p = poly(r);
   cvec r2 = roots(p);
   cout << "Roots, r = " << r << endl;
   cout << "Polynomial, p = " << p << endl;
   cout << "r = roots(p) = " << r2 << endl;

   r = randn_c(7);
   p = poly(r);
   r2 = roots(p);
   cout << "Roots, r = " << r << endl;
   cout << "Polynomial, p = " << p << endl;
   cout << "r = roots(p) = " << r2 << endl;

   cvec x = randn_c(9);
   cvec y = polyval(p, x);
   cout << "x = " << x << endl;
   cout << "polyval(p, x) = " << y << endl;
 }

 return 0;
}

#else

int main() { 
 cerr << "Error: LAPACK library is needed to run this test program" << endl; 
 return 1;
}

#endif
