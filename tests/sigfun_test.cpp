#include <iostream>
#include <itpp/itbase.h>

using std::cout;
using std::endl;

using namespace itpp;

int main()
{
  RNG_reset(12345);

  vec x = randn(10);
  vec y = randn(10);
  int max_lag;

  cout << "x = " << x << endl;
  cout << "y = " << y << endl;

  cout << "Testing cross correlation" << endl;
  cout << "--------------------------------------------------------------------------" << endl;

  cout << "xcorr(x,y) = " << xcorr(x,y) << endl;
  max_lag = 0; cout << "xcorr(x,y,"<< max_lag <<") = " << xcorr(x,y,max_lag) << endl;
  max_lag = 5; cout << "xcorr(x,y,"<< max_lag <<") = " << xcorr(x,y,max_lag) << endl;
  max_lag = -1; cout << "xcorr(x,y,"<< max_lag <<") = " << xcorr(x,y,max_lag) << endl;

  max_lag = 0; cout << "xcorr(x,y,"<< max_lag <<",biased) = " << xcorr(x,y,max_lag,"biased") << endl;
  max_lag = 5; cout << "xcorr(x,y,"<< max_lag <<",biased) = " << xcorr(x,y,max_lag,"biased") << endl;
  max_lag = -1; cout << "xcorr(x,y,"<< max_lag <<",biased) = " << xcorr(x,y,max_lag,"biased") << endl;

  max_lag = 0; cout << "xcorr(x,y,"<< max_lag <<",unbiased) = " << xcorr(x,y,max_lag,"unbiased") << endl;
  max_lag = 5; cout << "xcorr(x,y,"<< max_lag <<",unbiased) = " << xcorr(x,y,max_lag,"unbiased") << endl;
  max_lag = -1; cout << "xcorr(x,y,"<< max_lag <<",unbiased) = " << xcorr(x,y,max_lag,"unbiased") << endl;

  max_lag = 0; cout << "xcorr(x,y,"<< max_lag <<",coeff) = " << xcorr(x,y,max_lag,"coeff") << endl;
  max_lag = 5; cout << "xcorr(x,y,"<< max_lag <<",coeff) = " << xcorr(x,y,max_lag,"coeff") << endl;
  max_lag = -1; cout << "xcorr(x,y,"<< max_lag <<",coeff) = " << xcorr(x,y,max_lag,"coeff") << endl;

  cout << "Testing auto correlation" << endl;
  cout << "--------------------------------------------------------------------------" << endl;

  cout << "xcorr(x) = " << xcorr(x) << endl;
  max_lag = 0;  cout << "xcorr(x,"<< max_lag <<") = " << xcorr(x,max_lag) << endl;
  max_lag = 5;  cout << "xcorr(x,"<< max_lag <<") = " << xcorr(x,max_lag) << endl;
  max_lag = -1; cout << "xcorr(x,"<< max_lag <<") = " << xcorr(x,max_lag) << endl;

  max_lag = 0;  cout << "xcorr(x,"<< max_lag <<",biased) = " << xcorr(x,max_lag,"biased") << endl;
  max_lag = 5;  cout << "xcorr(x,"<< max_lag <<",biased) = " << xcorr(x,max_lag,"biased") << endl;
  max_lag = -1; cout << "xcorr(x,"<< max_lag <<",biased) = " << xcorr(x,max_lag,"biased") << endl;

  max_lag = 0;  cout << "xcorr(x,"<< max_lag <<",unbiased) = " << xcorr(x,max_lag,"unbiased") << endl;
  max_lag = 5;  cout << "xcorr(x,"<< max_lag <<",unbiased) = " << xcorr(x,max_lag,"unbiased") << endl;
  max_lag = -1; cout << "xcorr(x,"<< max_lag <<",unbiased) = " << xcorr(x,max_lag,"unbiased") << endl;

  max_lag = 0;  cout << "xcorr(x,"<< max_lag <<",coeff) = " << xcorr(x,max_lag,"coeff") << endl;
  max_lag = 5;  cout << "xcorr(x,"<< max_lag <<",coeff) = " << xcorr(x,max_lag,"coeff") << endl;
  max_lag = -1; cout << "xcorr(x,"<< max_lag <<",coeff) = " << xcorr(x,max_lag,"coeff") << endl;

  return 0;
}
