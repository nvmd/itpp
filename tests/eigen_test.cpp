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

int main(void)
{

  cout << "================================" << endl;
  cout << "  Test of eigenvalue routines" << endl;
  cout << "================================" << endl;

  {
    cout << "Real symmetric matrix" << endl;
    mat A = randn(5,5);
    A = transpose(A)*A; // make it symmetic
    mat V;
    vec d;
    eig_sym(A, d, V);

    cout << "A = " << A << endl;
    cout << "V = " << V << endl;
    cout << "d = " << d << endl;
    cout << "only d = " << eig_sym(A) << endl;
    cout << "norm(A*V-V*diag(d)) = " << norm(A*V-V*diag(d)) << endl;
  }

  {
    cout << endl << "Real non-symmetric matrix" << endl;
    mat A = randn(5,5);
    cmat V;
    cvec d;
    eig(A, d, V);

    cout << "A = " << A << endl;
    cout << "V = " << V << endl;
    cout << "d = " << d << endl;
    cout << "only d = " << eig(A) << endl;
    cout << "norm(A*V-V*diag(d)) = " << norm(A*V-V*diag(d)) << endl;
  }

  {
    cout << endl << "Complex hermitian matrix" << endl;
    cmat A = randn_c(5,5);
    A = transpose(conj(A))*A; // make it hermitian
    cmat V;
    vec d;
    eig_sym(A, d, V);

    cout << "A = " << A << endl;
    cout << "V = " << V << endl;
    cout << "d = " << d << endl;
    cout << "only d = " << eig_sym(A) << endl;
    cout << "norm(A*V-V*diag(d)) = " << norm(A*V-V*diag(d)) << endl;
  }

  {
    cout << endl << "Complex non-hermitian matrix" << endl;
    cmat A = randn_c(5,5);
    cmat V;
    cvec d;
    eig(A, d, V);

    cout << "A = " << A << endl;
    cout << "V = " << V << endl;
    cout << "d = " << d << endl;
    cout << "only d = " << eig(A) << endl;
    cout << "norm(A*V-V*diag(d)) = " << norm(A*V-V*diag(d)) << endl;
  }

  return 0;

}

#endif
