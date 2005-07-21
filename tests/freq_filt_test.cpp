/*
 * Test the Freq_Filt Class
 *
 * $Revision$
 *
 * $Date$
 */


#include <iostream>
#include <itpp/itbase.h>

using std::cout;
using std::endl;
using namespace itpp;

int main()
{
  vec b = "1 2 3 4";
  vec x(20);
  x.zeros();
  x(0) = 1;

  // Define a filter object for doubles
  Freq_Filt<double> FF(b,x.length());

  // Filter the data
  vec y = FF.filter(x);

  // Check the FFT and block sizes that were used
  int fftsize = FF.get_fft_size();
  int blksize = FF.get_blk_size();

  cout << fftsize << endl;
  cout << blksize << endl;

  cout << y << endl;

  // Test streaming mode
  x = linspace(0,10,100);
  Freq_Filt<double> FFS(b,x.length());
  vec y1 = FFS.filter(x(0,49),1);
  vec y2 = FFS.filter(x(50,99),1);

  cout << concat(y1,y2) << endl;

  return 0;

}
