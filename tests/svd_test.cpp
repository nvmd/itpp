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
int main() { cout << "LAPACK and CBLAS is needed for this test program" << endl; }
#else

int main(void)
{
  cout << "================================" << endl;
  cout << "      Test of svd routines" << endl;
  cout << "================================" << endl;
  {
    cout << "Real matrix" << endl;
    mat A = randn(5,5);
    mat U, V;
    vec S;
    svd(A, U, S, V);

    cout << "A = " << A << endl;
    cout << "U = " << U << endl;
    cout << "V = " << V << endl;
    cout << "S = " << S << endl;
    //cout << "A = U*diag(S)*V^T = " << U*diag(S)*transpose(V) << endl;
    cout << "only S = " << svd(A) << endl;

  }
  {
    cout << endl << "Complex matrix" << endl;
    cmat A = randn_c(5,5);
    cmat U, V;
    vec S;
    svd(A, U, S, V);

    cout << "A = " << A << endl;
    cout << "U = " << U << endl;
    cout << "V = " << V << endl;
    cout << "S = " << S << endl;
    //cout << "A = U*diag(S)*V^H = " << U*diag(S)*conj(transpose(V)) << endl;
    cout << "only S = " << svd(A) << endl;
  }

  return 0;

}

#endif
