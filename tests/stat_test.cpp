#include "itpp/itbase.h"


using std::cout;
using std::endl;
using namespace itpp;

int main()
{
  cout << "=================================" << endl;
  cout << "  Test of statistical routines" << endl;
  cout << "=================================" << endl;

  vec a = randn(5);

  cout << "a = " << a << endl;
  cout << "max(a) = " << max(a) << endl;
  cout << "min(a) = " << min(a) << endl;

  mat A = randn(5,5);
  cout << "A = " << A << endl << endl;

  cout << "max(A) = " << max(A) << endl;
  cout << "max(A,1) = " << max(A,1) << endl;
  cout << "max(A,2) = " << max(A,2) << endl;
  cout << "min(A) = " << min(A) << endl;
  cout << "min(A,1) = " << min(A,1) << endl;
  cout << "min(A,2) = " << min(A,2) << endl;

  cout << "norm(A) = " << norm(A) << endl;
  cout << "norm(A,2) = " << norm(A,2) << endl;
  cout << "norm(A,1) = " << norm(A,1) << endl;



  return 0;
}
