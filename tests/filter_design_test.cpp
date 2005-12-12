#include <itpp/itbase.h>

using std::cout;
using std::endl;
using std::complex;
using namespace itpp;

#ifndef HAVE_LAPACK
#define __THIS_PROGRAM_WILL_NOT_RUN__
#endif

#ifndef HAVE_CBLAS
#define __THIS_PROGRAM_WILL_NOT_RUN__
#endif

#ifdef __THIS_PROGRAM_WILL_NOT_RUN__

int main() 
{ 
  cout << "LAPACK and CBLAS is needed for this test program" << endl;
}

#else

int main()
{
  cout << "===================================" << endl;
  cout << "    Test of filter design routines" << endl;
  cout << "===================================" << endl;

  {
    cout << "Stabilisation of real filters" << endl;
    vec p = "0.7 3.0 -0.4";
    vec p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << p2 << endl;

    p = randn(7);
    p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << p2 << endl;

    p = poly(vec("1.1 0.7 0.2"));
    p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << p2 << endl;
  }

  {
    cout << "Stabilisation of complex filters" << endl;
    cvec p = "0.7 3.0 -0.4";
    cvec p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << p2 << endl;

    p = randn_c(7);
    p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << p2 << endl;

    p = poly(cvec("1.1 0.7 0.2"));
    p2 = polystab(p);
    cout << "Polynomial, p = " << p << endl;
    cout << "p2 = polystab(p) = " << p2 << endl;

    cvec a = randn_c(4);
    cvec b = randn_c(6);
    cvec h = freqz(b, a, 32);

    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << "freqz(b,a,32) = " << h << endl;
    
  }

  return 0;
}

#endif
