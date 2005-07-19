#include "itpp/itbase.h"

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
int main() { cout << "LAPACK and CBLAS is needed for this test program" << endl; }
#else

int main()
{

  cout << "================================" << endl;
  cout << "    Test of Cholesky routines" << endl;
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
      cout << "X = " << X << endl;
      if (!ok)
	cout << "matrix is not positive definite" << endl;
      else {
      cout << "F = " << F << endl;
      cout << "X = F^T * F = " << F.T() * F << endl;
      cout << "only F = " << chol(X) << endl;
      cout << "norm(e) = " << norm(X-F.T()*F) << endl;
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
      cout << "X = " << X << endl;

      if (!ok)
	cout << "matrix is not positive definite" << endl;
      else {
      cout << "F = " << F << endl;
      cout << "X = F^H * F = " << F.H() * F << endl;
      cout << "only F = " << chol(X) << endl;
      cout << "norm(e) = " << norm(X-F.H()*F) << endl;
      }
    }
  }

  return 0;

}

#endif
