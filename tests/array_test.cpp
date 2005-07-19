#include <iostream>
#include "itpp/base/array.h"

using namespace std;
using namespace itpp;

int main()
{
  Array<int> v1(10), v2(15);

  v1 = 1;
  v1(2) = 42;
  v2 = 5;

  cout << v1 << endl
       << v2 << endl;
	
  return 0;
}
