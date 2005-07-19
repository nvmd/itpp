#include "itpp/itbase.h"
#include "itpp/base/qr.h"

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
  cout << "=================================" << endl;
  cout << "Test of QR factorization routines" << endl;
  cout << "=================================" << endl;
  {
    cout << "QR of Real matrix" << endl;
    cout << "=================" << endl;
    mat A = randn(5,5);
    mat Q, R, e;

    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn(4,2);
    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;

    A = randn(2,4);
    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;
  }

  {
    cout << "QR of Real matrix with pivoting" << endl;
    cout << "===============================" << endl;
    mat A = randn(5,5);
    mat Q, R, e;
    bmat P;

    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn(4,2);
    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;

    A = randn(2,4);
    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;
  }

  {
    cout << "QR of Complex matrix" << endl;
    cout << "====================" << endl;
    cmat A = randn_c(5,5);
    cmat Q, R, e;

    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn_c(4,2);
    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;

    A = randn_c(2,4);
    qr(A, Q, R);
    e = A-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;
  }

  {
    cout << "QR of Complex matrix with pivoting" << endl;
    cout << "==================================" << endl;
    cmat A = randn_c(5,5);
    cmat Q, R, e;
    bmat P;

    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;

    // This does not give same sizes as matlab. Why???!!!!!
    A = randn_c(4,2);
    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;

    A = randn_c(2,4);
    qr(A, Q, R, P);
    e = A*to_mat(P)-Q*R;

    cout << "A = " << A << endl;
    cout << "Q = " << Q << endl;
    cout << "R = " << R << endl;
    cout << "P = " << P << endl << endl;
    cout << "norm(e) = " << norm(e) << endl << endl;
  }

  return 0;

}

#endif
