#include "itpp/itcomm.h"

using std::cout;
using std::endl;
using namespace itpp;

int main(void)
{
  GF a(8), b(8), c(8);

  a = 4;
  b = 2;

  c = a+b;

  cout << "a=" << a <<", b=" << b << endl;
  cout << "c=" << c << endl;

  return 0;

}
