#include <iostream>
#include <itpp/itbase.h>

using std::cout;
using std::endl;
using std::complex;
using namespace itpp;

int main()
{
   vec v1(10), v2(10), v3="1:4";
   cvec v4(4);

   v1 = 1.0;
   v2 = 5.0;
   v4(0) = complex<double>(1.2, 3.4);
   v4(1) = complex<double>(2.2, 1.6);
   v4(2) = complex<double>(3.2, 2.3);
   v4(3) = complex<double>(2.4, 5.4);

   cout << v1 << " + " << v2 << " = " << v1+v2 << endl;

   cout << "norm(" << v3 << ") = " << norm(v3) << endl;
   cout << "norm(" << v3 << ",3) = " << norm(v3,3) << endl;
   cout << "energy(" << v3 << ") = " << energy(v3) << endl;
   cout << "variance(" << v3 << ") = " << variance(v3) << endl;

   cout << "norm(" << v4 << ") = " << norm(v4) << endl;
   cout << "norm(" << v4 << ",3) = " << norm(v4,3) << endl;
   cout << "energy(" << v4 << ") = " << energy(v4) << endl;
   cout << "variance(" << v4 << ") = " << variance(v4) << endl;

   cout << "hamming(4)" << hamming(4) << endl
	 << "hamming(5)" << hamming(5) << endl
	 << "triang(4)" << triang(4) << endl
	 << "triang(5)" << triang(5) << endl
	 << "sqrt_win(4)" << sqrt_win(4) << endl
	 << "sqrt_win(5)" << sqrt_win(5) << endl;

   v1 = ones(5);
   cout << "v1 = " << v1 << endl;
   v1.shift_left(2.0,2);
   cout << "v1.shift_left(2.0,2): " << v1 << endl;
   v1.shift_right(3.0);
   cout << "v1.shift_left(3.0): " << v1 << endl;
   v1.shift_left(vec("4 5 6"));
   cout << "v1.shift_left(vec(\"4 5 6\")): " << v1 << endl;
   v1.shift_right(vec("7 8"));
   cout << "v1.shift_right(vec(\"7 8\")): " << v1 << endl;

   //Test of rem:
   vec x = "1.0 2.0 3.4 -4.5 6.7";
   double y = 0.76;
   cout << "x = " << x << endl;
   cout << "y = " << y << endl;
   cout << "rem(x,y) = " << rem(x,y) << endl;
   cout << "rem(10,x) = " << rem(10,x) << endl;
   mat xx = "1.0 2.3;4.5 -6.7";
   cout << "xx = " << xx << endl;
   cout << "rem(xx,y) = " << rem(xx,y) << endl;
   cout << "rem(10,xx) = " << rem(10,xx) << endl;

   //Test of to_str:
   int    a = 1;               cout << to_str(a) << endl;
   long   b = 2;               cout << to_str(b) << endl;
   short  c = 3;               cout << to_str(c) << endl;
   double d = 4.0;             cout << to_str(d) << endl;
   ivec   e = "5 6 7";         cout << to_str(e) << endl;
   complex<double> f(8.0,9.0); cout << to_str(f) << endl;

   vec aa = "1.0   ,  2.0 3 -4.5 6.7";
   cout << "aa = " << aa << endl;

   bvec bb = "1 0 ,  1 0";
   cout << "bb = " << bb << endl;

   mat A1, A2, A3, A4;
   A1 = "1.0   ,  2.0;3.0 4.0";   cout << "A1 = " << A1 << endl;
   A2 = "1.0 2.0; 3.0 4.0";  cout << "A2 = " << A2 << endl;
   A3 = "1.0 2.0    ;   3.0 4.0";  cout << "A3 = " << A3 << endl;
   A4 = "1.0 2.0 ; 3.0 4.0"; cout << "A4 = " << A4 << endl;

   return 0;
}
