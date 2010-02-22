#include <itpp/itbase.h>

using namespace std;
using namespace itpp;

int main(void)
{

  cout << "================================" << endl;
  cout << "    Test of bessel functions " << endl;
  cout << "================================" << endl;

  vec x = linspace(0.01, 10, 20);

  cout << "x = " << x << endl;

  cout << "besselj(0, x) = " << fixed << besselj(0, x) << endl;
  cout << "besselj(1, x) = " << fixed << besselj(1, x) << endl;
  cout << "besselj(5, x) = " << fixed << besselj(5, x) << endl;
  cout << "besselj(0.3, x) = " << fixed << besselj(0.3, x) << endl;
  cout << "besselj(1.7, x) = " << fixed << besselj(1.7, x) << endl;
  cout << "besselj(5.3, x) = " << fixed << besselj(5.3, x) << endl;

  cout << "bessely(0, x) = " << fixed << bessely(0, x) << endl;
  cout << "bessely(1, x) = " << fixed << bessely(1, x) << endl;
  cout << "bessely(5, x) = " << fixed << bessely(5, x) << endl;
  cout << "bessely(0.3, x) = " << fixed << bessely(0.3, x) << endl;
  cout << "bessely(1.7, x) = " << fixed << bessely(1.7, x) << endl;
  cout << "bessely(5.3, x) = " << fixed << bessely(5.3, x) << endl;

  cout << "besseli(0, x) = " << fixed << besseli(0, x) << endl;
  cout << "besseli(1, x) = " << fixed << besseli(1, x) << endl;
  cout << "besseli(5, x) = " << fixed << besseli(5, x) << endl;
  cout << "besseli(0.3, x) = " << fixed << besseli(0.3, x) << endl;
  cout << "besseli(1.7, x) = " << fixed << besseli(1.7, x) << endl;
  cout << "besseli(5.3, x) = " << fixed << besseli(5.3, x) << endl;

  cout << "besselk(0, x) = " << fixed << besselk(0, x) << endl;
  cout << "besselk(1, x) = " << fixed << besselk(1, x) << endl;
  cout << "besselk(5, x) = " << fixed << besselk(5, x) << endl;

  return 0;

}
