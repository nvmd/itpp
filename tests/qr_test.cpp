/*!
 * \file 
 * \brief QR factorisation test program
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


#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

int main()
{
  cout << "=================================" << endl;
  cout << "Test of QR factorization routines" << endl;
  cout << "=================================" << endl;

  {
    cout << "QR of Real matrix" << endl;
    cout << "=================" << endl;
    mat A = randn(5,5);
    mat Q, R, e;

    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn(4,2);
    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;

    A = randn(2,4);
    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;
  }

  {
    cout << "QR of Real matrix with pivoting" << endl;
    cout << "===============================" << endl;
    mat A = randn(5,5);
    mat Q, R, e;
    bmat P;

    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn(4,2);
    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;

    A = randn(2,4);
    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;
  }

  {
    cout << "QR of Complex matrix" << endl;
    cout << "====================" << endl;
    cmat A = randn_c(5,5);
    cmat Q, R, e;

    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn_c(4,2);
    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;

    A = randn_c(2,4);
    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;
  }

  {
    cout << "QR of Complex matrix with pivoting" << endl;
    cout << "==================================" << endl;
    cmat A = randn_c(5,5);
    cmat Q, R, e;
    bmat P;

    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn_c(4,2);
    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;

    A = randn_c(2,4);
    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << round_to_zero(A) << endl;
    cout << "Q = " << round_to_zero(Q) << endl;
    cout << "R = " << round_to_zero(R) << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << round_to_zero(norm(e), 1e-14) << endl << endl;
  }

  return 0;

}

#else

int main() { 
  cerr << "Error: LAPACK (or MKL) is needed for this test program" << endl;
  return 1;
}

#endif
