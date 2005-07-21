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
int main() { cout << "LAPACK and CBLAS is needed for this test program" << endl; }
#else

int main()
{

  cout << "===================================" << endl;
  cout << "    Test of Determinant routines" << endl;
  cout << "===================================" << endl;

  {
    cout << "Real matrix" << endl;
    mat X = randn(5,5);
    double d;
    d = det(X);
    cout << "X = " << X << endl;
    cout << "det(X) = " << d << endl;

    X = randn(5,5);
    d = det(X);
    cout << "X = " << X << endl;
    cout << "det(X) = " << d << endl;
  }

  {
    cout << endl << "Complex matrix" << endl;
    cmat X = randn_c(5,5);
    complex<double> d;
    d = det(X);
    cout << "X = " << X << endl;
    cout << "det(X) = " << d << endl;

    X = randn_c(5,5);
    d = det(X);
    cout << "X = " << X << endl;
    cout << "det(X) = " << d << endl;
  }

  return 0;

}

#endif
