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
  cout << "=================================" << endl;
  cout << "Test of LU factorization routines" << endl;
  cout << "=================================" << endl;
  {
    cout << "Real matrix" << endl;
    mat X = randn(5,5);
    mat L, U;
    ivec p;
    lu(X, L, U, p);


    mat P=to_mat(permutation_matrix(p));
    cout << "X = " << X << endl;
    cout << "L = " << L << endl;
    cout << "U = " << U << endl;
    cout << "p = " << p << endl;
    cout << "P = " << P << endl;
    cout << "X = P^T * L * U = " << transpose(P)*L*U << endl;
  }
  {
    cout << "Complex matrix" << endl;
    cmat X = randn_c(5,5);
    cmat L, U;
    ivec p;
    lu(X, L, U, p);


    mat P=to_mat(permutation_matrix(p));
    cout << "X = " << X << endl;
    cout << "L = " << L << endl;
    cout << "U = " << U << endl;
    cout << "p = " << p << endl;
    cout << "P = " << P << endl;
    cout << "X = P^T * L * U = " << transpose(P)*L*U << endl;
  }

  return 0;

}

#endif
