#include <itpp/itbase.h>


using std::cout;
using std::endl;
using namespace itpp;

int main()
{
  cout << "=================================" << endl;
  cout << "    Test of matfunc routines" << endl;
  cout << "=================================" << endl;

  vec a = randn(5);
  cout << "a = " << a << endl;
  cout << "sum(a) = " << sum(a) << endl;
  cout << "cumsum(a) = " << cumsum(a) << endl;
  cout << "prod(a) = " << prod(a) << endl;
  cout << "sum_sqr(a) = " << sum_sqr(a) << endl << endl;

  mat A = randn(5,5);
  cout << "A = " << A << endl << endl;

  cout << "sum(A) = " << sum(A) << endl;
  cout << "sum(A,1) = " << sum(A,1) << endl;
  cout << "sum(A,2) = " << sum(A,2) << endl << endl;

  cout << "cumsum(A) = " << cumsum(A) << endl;
  cout << "cumsum(A,1) = " << cumsum(A,1) << endl;
  cout << "cumsum(A,2) = " << cumsum(A,2) << endl << endl;

  cout << "prod(A) = " << prod(A) << endl;
  cout << "prod(A,1) = " << prod(A,1) << endl;
  cout << "prod(A,2) = " << prod(A,2) << endl << endl;

  cout << "sum_sqr(A) = " << sum_sqr(A) << endl;
  cout << "sum_sqr(A,1) = " << sum_sqr(A,1) << endl;
  cout << "sum_sqr(A,2) = " << sum_sqr(A,2) << endl;

  return 0;
}
