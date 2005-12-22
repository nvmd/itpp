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

int main()
{
  cout << "==========================================" << endl;
  cout << "   Solving linear systems of equations" << endl;
  cout << "==========================================" << endl;

  {
    cout << "Real systems:" << endl << endl;
    mat A, B, X;
    vec b, x;

    A = randn(4,4);
    b = randn(4);
    x = ls_solve(A, b);

    cout << "Square system: Ax=b" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "b=" << b << endl
	 << "x=" << x << endl << endl;

    A = randn(4,4);
    B = randn(4,2);
    X = ls_solve(A, B);

    cout << "Square system: AX=B" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "B=" << B << endl
	 << "X=" << X << endl << endl;



    A = randn(4,4); A = A.transpose()*A;
    b = randn(4);
    x = ls_solve(A, b);

    cout << "Square system (chol): Ax=b" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "b=" << b << endl
	 << "x=" << x << endl << endl;


    A = randn(4,4); A = A.transpose()*A;
    B = randn(4,2);
    X = ls_solve(A, B);

    cout << "Square system (Chol): AX=B" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "B=" << B << endl
	 << "X=" << X << endl << endl;



    A = randn(4,2);
    b = randn(4);
    x = ls_solve_od(A, b);

    cout << "Overdetermined system: Ax=b" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "b=" << b << endl
	 << "x=" << x << endl << endl;

    A = randn(4,2);
    B = randn(4,3);
    X = ls_solve_od(A, B);

    cout << "Overdetermined system: AX=B" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "B=" << B << endl
	 << "X=" << X << endl << endl;

    A = randn(2,4);
    b = randn(2);
    x = ls_solve_ud(A, b);

    cout << "Underdetermined system: Ax=b" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "b=" << b << endl
	 << "x=" << x << endl << endl;

    A = randn(2,4);
    B = randn(2,3);
    X = ls_solve_ud(A, B);

    cout << "Underdetermined system: AX=B" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "B=" << B << endl
	 << "X=" << X << endl << endl;

  }


  {
    cout << "Complex systems:" << endl << endl;
    cmat A, B, X;
    cvec b, x;

    A = randn_c(4,4);
    b = randn_c(4);
    x = ls_solve(A, b);

    cout << "Square system: Ax=b" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "b=" << b << endl
	 << "x=" << x << endl << endl;

    A = randn_c(4,4);
    B = randn_c(4,2);
    X = ls_solve(A, B);

    cout << "Square system: AX=B" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "B=" << B << endl
	 << "X=" << X << endl << endl;



    A = randn_c(4,4); A = A.transpose()*A;
    b = randn_c(4);
    x = ls_solve(A, b);

    cout << "Square system (chol): Ax=b" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "b=" << b << endl
	 << "x=" << x << endl << endl;


    A = randn_c(4,4); A = A.transpose()*A;
    B = randn_c(4,2);
    X = ls_solve(A, B);

    cout << "Square system (Chol): AX=B" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "B=" << B << endl
	 << "X=" << X << endl << endl;



    A = randn_c(4,2);
    b = randn_c(4);
    x = ls_solve_od(A, b);

    cout << "Overdetermined system: Ax=b" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "b=" << b << endl
	 << "x=" << x << endl << endl;

    A = randn_c(4,2);
    B = randn_c(4,3);
    X = ls_solve_od(A, B);

    cout << "Overdetermined system: AX=B" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "B=" << B << endl
	 << "X=" << X << endl << endl;

    A = randn_c(2,4);
    b = randn_c(2);
    x = ls_solve_ud(A, b);

    cout << "Underdetermined system: Ax=b" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "b=" << b << endl
	 << "x=" << x << endl << endl;

    A = randn_c(2,4);
    B = randn_c(2,3);
    X = ls_solve_ud(A, B);

    cout << "Underdetermined system: AX=B" << endl
	 << "============================" << endl
	 << "A=" << A << endl
	 << "B=" << B << endl
	 << "X=" << X << endl << endl;

  }
  return 0;

}

#endif
