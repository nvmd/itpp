/*!
 * \file 
 * \brief Cholesky factorisation test program
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


#if (defined(HAVE_LAPACK) && defined(HAVE_CBLAS)) || defined(HAVE_MKL)

int main()
{
  cout << "================================" << endl;
  cout << "    Test of Cholesky routines   " << endl;
  cout << "================================" << endl;

  {
    cout << "Real matrix" << endl;
    mat X, F;
    bool ok = false;

    while (!ok) {
      X = randn(5,5);
      X = X.T()*X; // create a symmetric matrix
      // Make diagonal real and positive
      for (int i=0; i<X.cols(); i++)
	X(i,i) = std::abs(X(i,i));

      ok = chol(X, F);      
      cout << "X = " << round_to_zero(X) << endl;
      if (!ok)
	cout << "matrix is not positive definite" << endl;
      else {
	cout << "F = " << round_to_zero(F) << endl;
	cout << "X = F^T * F = " << round_to_zero(F.T() * F) << endl;
	cout << "only F = " << round_to_zero(chol(X)) << endl;
	cout << "norm(e) = " << round_to_zero(norm(X-F.T()*F)) << endl;
      }
    }
  }

  {
    cout << endl << "Complex matrix" << endl;
    cmat X, F;
    bool ok = false;

    while (!ok) {
      X = randn_c(5,5);
      X = X.H()*X; // create a symmetric matrix
      // Make diagonal real and positive
      for (int i=0; i<X.cols(); i++)
	X(i,i) = std::abs(real(X(i,i)))
;
      ok = chol(X, F);      
      cout << "X = " << round_to_zero(X) << endl;

      if (!ok)
	cout << "matrix is not positive definite" << endl;
      else {
	cout << "F = " << round_to_zero(F) << endl;
	cout << "X = F^H * F = " << round_to_zero(F.H() * F) << endl;
	cout << "only F = " << round_to_zero(chol(X)) << endl;
	cout << "norm(e) = " << round_to_zero(norm(X-F.H()*F)) << endl;
      }
    }
  }

  return 0;
}

#else

int main() { 
  cerr << "Error: LAPACK and CBLAS (or MKL) are needed for this test program" << endl;
  return 1;
}

#endif
