/*!
 * \file 
 * \brief SVD decomposition test program
 * \author Tony Ottosson and Adam Piatyszek
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

#include <itpp/itstat.h>

using namespace itpp;
using namespace std;


#if defined(HAVE_LAPACK)

int main()
{
  cout << "================================" << endl;
  cout << "      Test of svd routines      " << endl;
  cout << "================================" << endl;

  {
    cout << "Real matrix" << endl;
    mat A = randn(5, 5);
    mat U, V;
    vec S;
    svd(A, U, S, V);

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - U*diag(S)*V^T) = " 
	 << round_to_zero(norm(A - U * diag(S) * transpose(V))) << endl;

  }
  {
    cout << endl << "Complex matrix" << endl;
    cmat A = randn_c(5, 5);
    cmat U, V;
    vec S;
    svd(A, U, S, V);

    cout << "A = " << round_to_zero(A) << endl;
    cout << "norm(A - U*diag(S)*V^H) = " 
	 << round_to_zero(norm(A - U * diag(S) * hermitian_transpose(V))) 
	 << endl;
  }

  return 0;
}

#else

int main() { 
  cerr << "Error: LAPACK library is needed to run this test program" << endl;
  return 1;
}

#endif
