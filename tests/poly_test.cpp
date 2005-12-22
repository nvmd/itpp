#include <itpp/itbase.h>

using std::cout;
using std::endl;
using std::complex;
using namespace itpp;

#ifdef NO_LAPACK
#define __THIS_PROGRAM_WILL_NOT_RUN__
#endif

#ifdef NO_CBLAS
#define __THIS_PROGRAM_WILL_NOT_RUN__
#endif

#ifdef __THIS_PROGRAM_WILL_NOT_RUN__
int main() { cout << "LAPACK and CBLAS (or MKL) are needed for this test program" << endl; }
#else

int main()
{

  cout << "===================================" << endl;
  cout << "    Test of polynomial routines" << endl;
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

#endif
