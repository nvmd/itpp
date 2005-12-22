#include <itpp/itbase.h>

using std::cout;
using std::endl;
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

int main(void)
{
  cout << "======================================" << endl;
  cout << "    Test of Matrix inversion routines" << endl;
  cout << "======================================" << endl;
  {
    cout << "Real matrix" << endl;

    mat X = randn(5,5), Y;
    Y = inv(X);
    cout << "X = " << X << endl;
    cout << "inv(X) = " << Y << endl;

    X = randn(5,5);
    Y = inv(X);
    cout << "X = " << X << endl;
    cout << "inv(X) = " << Y << endl;
  }
  {
    cout << endl << "Complex matrix" << endl;

    cmat X = randn_c(5,5), Y;
    Y = inv(X);
    cout << "X = " << X << endl;
    cout << "inv(X) = " << Y << endl;

    X = randn_c(5,5);
    Y = inv(X);
    cout << "X = " << X << endl;
    cout << "inv(X) = " << Y << endl;
  }

  return 0;

}

#endif
