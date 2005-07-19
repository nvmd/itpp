#include "itpp/itbase.h"

using std::cout;
using std::endl;
using namespace itpp;

int main()
{
  RNG_reset(4357U);
  cout << "randn(5) = " << randn(5) << endl;

  return 0;
}
