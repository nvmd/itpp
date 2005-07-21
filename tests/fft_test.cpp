#include <itpp/itbase.h>

using std::cout;
using std::endl;
using namespace itpp;

int main()
{
  int N=pow2i(16), NN=1000, i;
  vec a(N), b(N);
  double *p1=a._data(),*p2=b._data();

  //Remove comments below to measure execution times.
  //Real_Timer time;

  //time.tic();
  for (i=0; i<NN; i++) {
    memcpy(p1,p2,4*N);
  }
  //time.toc();

  N=17, NN=100, i;

  vec in(N), in_final;
  cvec out, out_real, in_complex, inv_out;

  in = randn(N);

  cout << "fft:" << endl;
  //time.tic();
  for (i=0; i<NN; i++) {
    in_complex = to_cvec(in);
    fft(in_complex, out);
  }
  //time.toc();

  cout << "Real fft:" << endl;
  //time.tic();
  for (i=0; i<NN; i++) {
    fft_real(in, out_real);
  }
  //time.toc();

  cout << "Real ifft:" << endl;
  //time.tic();
  for (i=0; i<NN; i++) {
    ifft_real(out_real, in_final);
  }
  //time.toc();

  cout << "Old ifft:" << endl;
  //time.tic();
  for (i=0; i<NN; i++) {
    inv_out = ifft(out);
  }
  //time.toc();

  return 0;

}
