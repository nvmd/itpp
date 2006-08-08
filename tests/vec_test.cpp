/*!
 * \file 
 * \brief Vector class test program
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
 * Copyright (C) 1995-2006  (see AUTHORS file for a list of contributors)
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

int main()
{
  int N = 33;
  vec a = randn(N), b = randn(N), r, v;
  ivec iv;
  svec sv;
  cvec cv;
  bvec bv;
  double c = randn(), s , y;
  mat M;

  cout << "a = " << a << endl;
  cout << "b = " << b << endl;
  cout << "c = " << c << endl;

  r = a+b;
  cout << "r = a+b = " << r << endl;

  r = a+c;
  cout << "r = a+c = " << r << endl;

  r = c+a;
  cout << "r = c+a = " << r << endl;

  r = a-b;
  cout << "r = a-b = " << r << endl;

  r = a-c;
  cout << "r = a-c = " << r << endl;

  r = c-a;
  cout << "r = c-a = " << r << endl;

  r = -a;
  cout << "r = -a = " << r << endl;

  s = a*b;
  cout << "s = a*b = " << s << endl;

  s = dot(a,b);
  cout << "s = dot(a,b) = " << s << endl;

  M = outer_product(a,b);
  cout << "M = outer_product(a,b) = " << M << endl;

  r = a*c;
  cout << "r = a*c = " << r << endl;

  r = c*a;
  cout << "r = c*a = " << r << endl;

  r = elem_mult(a,b);
  cout << "r = elem_mult(a,b) = " << r << endl;

  r = a/c;
  cout << "r = a/c = " << r << endl;

  r = elem_div(a,b);
  cout << "r = elem_div(a,b) = " << r << endl;

  r = concat(a,c);
  cout << "r = concat(a,c) = " << r << endl;

  r = concat(c,a);
  cout << "r = concat(c,a) = " << r << endl;

  r = concat(a,b);
  cout << "r = concat(a,b) = " << r << endl;

  r = concat(a,b,a);
  cout << "r = concat(a,b,a) = " << r << endl;

  s = a.size();
  cout << "Testing Vec<T>.size(): s = " << s << endl;

  r.set_length(17,false);
  cout << "Testing Vec<T>.set_length(): r.size() = " << r.size() << endl;

  r.zeros();
  cout << "Testing Vec<T>.zeros(): r = " << r << endl;

  r.ones();
  cout << "Testing Vec<T>.ones(): r = " << r << endl << endl;

  // Test vectror initialisation with string
  v = "23.3 1232.7 0.111 1.525 0.333";
  cout << "Testing double vector initialisation with: \"23.3 1232.7 0.111 1.525 0.333\":" 
       << endl << "v = " << v << endl;

  v = "-10.000 :.5:-4.5";
  cout << "Testing double vector initialisation with: \"-10.000 :.5:-4.5\":" 
       << endl << "v = " << v << endl;

  iv = "0xA : -010";
  cout << "Testing int vector initialisation with: \"0xA : -010\":" 
       << endl << "iv = " << iv << endl;

  sv = "3 0xF -10, 0133 0177, 0x0 ";
  cout << "Testing short int vector initialisation with: \"3 0xF -10, 0133 0177, 0x0 \":" 
       << endl << "sv = " << sv << endl;

  cv = " (0.3, 0.4)  .2-.01i, 1e-3+0.25i";
  cout << "Testing complex vector initialisation with: \" (0.3, 0.4)  .2-.01i, 1e-3+0.25i\":" 
       << endl << "cv = " << cv << endl;

  bv = "1 1 0,1  1  ,  0 ,1  ";
  cout << "Testing bit vector initialisation with: \"1 1 0,1  1  ,  0 ,1  \":" 
       << endl << "bv = " << bv << endl << endl;

  // Test of shifts
  v = ones(5);
  cout << "v = " << v << endl;
  v.shift_left(2.0,2);
  cout << "v.shift_left(2.0,2): " << v << endl;
  v.shift_right(3.0);
  cout << "v.shift_right(3.0): " << v << endl;
  v.shift_left(vec("4 5 6"));
  cout << "v.shift_left(vec(\"4 5 6\")): " << v << endl;
  v.shift_right(vec("7 8"));
  cout << "v.shift_right(vec(\"7 8\")): " << v << endl << endl;

  // Test of rem:
  v = "1.0 2.0 3.4 -4.5 6.7";
  y = 0.76;
  cout << "v = " << v << endl;
  cout << "y = " << y << endl;
  cout << "rem(v, y) = " << rem(v, y) << endl;
  cout << "rem(10, v) = " << rem(10, v) << endl;
  M = "1.0 2.3; 4.5 -6.7";
  cout << "M = " << M << endl;
  cout << "rem(M, y) = " << rem(M, y) << endl;
  cout << "rem(10, M) = " << rem(10, M) << endl << endl;

  // Test of all and any:
  bvec b1 = "0 0 0 0 0 0 0 1 0 0";
  bvec b2 = "0 0 0 0 0 0 0 0 0 0";
  bvec b3 = "1 1 1 1 1 1 1 1 1 1 1 1 1";
  bvec b4 = "1 1 1 1 1 1 1 1 1 1 1 0 1";

  cout << "any(b1) = " << any(b1) << endl;
  cout << "any(b2) = " << any(b2) << endl;
  cout << "all(b3) = " << all(b3) << endl;
  cout << "all(b4) = " << all(b4) << endl;

  return 0;
}
