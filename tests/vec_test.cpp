//
// vec_test.cpp
//
// $Id$
//

/* 
   Test elementary operations within the Vec class using type double
   The program vec_test_file.cpp produces the test file "vec_test.it".
   Created by Tony Ottosson, 2001-10-11.
*/

#include <iostream>
#include "itpp/itbase.h"

using std::cout;
using std::endl;
using namespace itpp;

int main()
{
  int N = 33;
  vec a = randn(N), b = randn(N), r;
  double c = randn(), s;
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

  r="23.3 1232.7 0.111 1.525 0.333";
  cout << "Testing to set vector by string: r = " << r << endl;

  s = a.size();
  cout << "Testing Vec<T>.size(): s = " << s << endl;

  r.set_length(17,false);
  cout << "Testing Vec<T>.set_length(): r.size() = " << r.size() << endl;

  r.zeros();
  cout << "Testing Vec<T>.zeros(): r = " << r << endl;

  r.ones();
  cout << "Testing Vec<T>.ones(): r = " << r << endl;

  //Test of all and any:
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
