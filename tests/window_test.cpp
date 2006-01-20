#include <itpp/itbase.h>

using std::cout;
using std::endl;
using namespace itpp;

int main(void)
{

  cout << "================================" << endl;
  cout << "    Test of window functions " << endl;
  cout << "================================" << endl;

  cout << "hamming(32) = " << hamming(32) << endl;
  cout << "hamming(128) = " << hamming(128) << endl;

  cout << "hanning(32) = " << hanning(32) << endl;
  cout << "hanning(128) = " << hanning(128) << endl;

  cout << "hann(32) = " << hann(32) << endl;
  cout << "hann(128) = " << hann(128) << endl;

  cout << "blackman(32) = " << blackman(32) << endl;
  cout << "blackman(128) = " << blackman(128) << endl;

  cout << "triang(32) = " << triang(32) << endl;
  cout << "triang(128) = " << triang(128) << endl;

  return 0;

}
