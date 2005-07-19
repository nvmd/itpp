#include "itpp/itbase.h"

using std::cout;
using std::endl;
using namespace itpp;

int main(void)
{

  cout << "================================" << endl;
  cout << "    Test of bessel functions " << endl;
  cout << "================================" << endl;

  vec x = linspace(0.01, 10, 20);

  cout << "x = " << x << endl;

  cout << "besselj(0, x) = " << besselj(0, x) << endl;
  cout << "besselj(1, x) = " << besselj(1, x) << endl;
  cout << "besselj(5, x) = " << besselj(5, x) << endl;
  cout << "besselj(0.3, x) = " << besselj(0.3, x) << endl;
  cout << "besselj(1.7, x) = " << besselj(1.7, x) << endl;
  cout << "besselj(5.3, x) = " << besselj(5.3, x) << endl;

  cout << "bessely(0, x) = " << bessely(0, x) << endl;
  cout << "bessely(1, x) = " << bessely(1, x) << endl;
  cout << "bessely(5, x) = " << bessely(5, x) << endl;
  cout << "bessely(0.3, x) = " << bessely(0.3, x) << endl;
  cout << "bessely(1.7, x) = " << bessely(1.7, x) << endl;
  cout << "bessely(5.3, x) = " << bessely(5.3, x) << endl;

  cout << "besseli(0, x) = " << besseli(0, x) << endl;
  cout << "besseli(1, x) = " << besseli(1, x) << endl;
  cout << "besseli(5, x) = " << besseli(5, x) << endl;
  cout << "besseli(0.3, x) = " << besseli(0.3, x) << endl;
  cout << "besseli(1.7, x) = " << besseli(1.7, x) << endl;
  cout << "besseli(5.3, x) = " << besseli(5.3, x) << endl;

  cout << "besselk(0, x) = " << besselk(0, x) << endl;
  cout << "besselk(1, x) = " << besselk(1, x) << endl;
  cout << "besselk(5, x) = " << besselk(5, x) << endl;

  return 0;

}
