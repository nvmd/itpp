#include <itpp/itbase.h>

using std::cout;
using std::endl;
using namespace itpp;


#ifdef NO_FFTW
int main() { cout << "FFTW is needed for this test program" << endl; }
#else

int main()
{
  cout << "=============================================================" << endl;
  cout << "                    Test of Transforms                       " << endl;
  cout << "=============================================================" << endl;

  {
    cout << "Real vector: fft_real(x,y), ifft(y,z)" << endl;
    int N=16;

    vec x, z;
    cvec y;
    
    x = randn(N);
    fft_real(x, y);
    ifft_real(y, z);
    
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
  }
  {
    cout << "Complex vector: fft_real(x,y), ifft(y,z)" << endl;
    int N=16;

    cvec x, y, z;
    
    x = randn_c(N);
    fft(x, y);
    ifft(y, z);
    
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;

  }
}


#endif




