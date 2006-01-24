/*!
 * \file 
 * \brief Test program of various functions on vectors and matrices
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

using namespace std;
using namespace itpp;

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

int main()
{
  cout << "=================================" << endl;
  cout << "    Test of matfunc routines     " << endl;
  cout << "=================================" << endl;

  vec a = randn(5);
  cout << "a = " << a << endl;
  cout << "sum(a) = " << sum(a) << endl;
  cout << "cumsum(a) = " << cumsum(a) << endl;
  cout << "prod(a) = " << prod(a) << endl;
  cout << "sum_sqr(a) = " << sum_sqr(a) << endl << endl;

  mat A = randn(5,5);
  cout << "A = " << A << endl << endl;

  cout << "sum(A) = " << sum(A) << endl;
  cout << "sum(A,1) = " << sum(A,1) << endl;
  cout << "sum(A,2) = " << sum(A,2) << endl << endl;

  cout << "cumsum(A) = " << cumsum(A) << endl;
  cout << "cumsum(A,1) = " << cumsum(A,1) << endl;
  cout << "cumsum(A,2) = " << cumsum(A,2) << endl << endl;

  cout << "prod(A) = " << prod(A) << endl;
  cout << "prod(A,1) = " << prod(A,1) << endl;
  cout << "prod(A,2) = " << prod(A,2) << endl << endl;

  cout << "sum_sqr(A) = " << sum_sqr(A) << endl;
  cout << "sum_sqr(A,1) = " << sum_sqr(A,1) << endl;
  cout << "sum_sqr(A,2) = " << sum_sqr(A,2) << endl << endl;

  cout << "Kronecker test" << endl;
  mat X = to_mat(randi(2, 2, 1, 4));
  mat Y = randn(3, 3);
  cout << "X = " << X << endl;
  cout << "Y = " << Y << endl;
  cout << "kron(X, Y) = " << kron(X, Y) << endl << endl;

  cout << "sqrtm of a real matrix" << endl;
  A = randn(3,3);
  cmat A_sqrtm = sqrtm(A);
  cout << "A = " << A << endl;
  cout << "norm(sqrtm(A) * sqrtm(A) - A) = " 
       << round_to_zero(norm(A_sqrtm * A_sqrtm - to_cmat(A)), 1e-14) 
       << endl << endl;

  cout << "sqrtm of a complex matrix" << endl;
  cmat B = randn_c(3,3);
  cmat B_sqrtm = sqrtm(B);
  cout << "B = " << B << endl;
  cout << "norm(sqrtm(B) * sqrtm(B) - B) = " 
       << round_to_zero(norm(B_sqrtm * B_sqrtm - B), 1e-14) << endl;

  return 0;
}

#else

int main() { 
  cerr << "Error: LAPACK (or MKL) is needed for this test program" << endl; 
  return 1;
}

#endif
