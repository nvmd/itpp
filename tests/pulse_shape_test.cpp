#include <iomanip>
#include <itpp/itcomm.h>

using namespace std;
using namespace itpp;

int main()
{
  cout.setf(ios::fixed);
  cout << "======================================" << endl;
  cout << "    Test of pulse shaping routines    " << endl;
  cout << "======================================" << endl << endl;

  BPSK bpsk;
  vec symbols, samples, rsymbols;
  Root_Raised_Cosine<double> rrc_tx(0.5), rrc_rx(0.5);
  Raised_Cosine<double> rc_tx(0.5);

  bpsk.modulate_bits(randb(20), symbols);
  samples = rrc_tx.shape_symbols(symbols);
  rsymbols = rrc_rx.shape_samples(samples);

 	cout << "*** Root Raised Cosine; real input ***" << endl << endl;
  cout << "pulse, RRC = " << setprecision(6) << rrc_tx.get_pulse_shape() 
    << endl << endl;
  cout << "symbols = " << setprecision(0) << symbols << endl << endl;
  cout << "samples = " << setprecision(6) << samples << endl << endl;
  cout << "received symbols =" << setprecision(6) << rsymbols << endl << endl;

  samples = rc_tx.shape_symbols(symbols);

  cout << "*** Raised Cosine; real input ***" << endl << endl;
  cout << "pulse, RC = " << setprecision(6) << rc_tx.get_pulse_shape() 
    << endl << endl;
  cout << "symbols = " << setprecision(0) << symbols << endl << endl;
  cout << "samples = " << setprecision(6) << samples << endl << endl;

  QPSK qpsk;
  cvec csymbols, csamples, crsymbols;
  Root_Raised_Cosine<complex<double> > crrc_tx(0.5), crrc_rx(0.5);
  Raised_Cosine<complex<double> > crc_tx(0.5);

  qpsk.modulate_bits(randb(40), csymbols);
  csamples = crrc_tx.shape_symbols(csymbols);
  crsymbols = crrc_rx.shape_samples(csamples);

  cout << "*** Root Raised Cosine; complex input ***" << endl << endl;
  cout << "pulse, RRC = " << setprecision(6) << crrc_tx.get_pulse_shape() 
    << endl << endl;
  cout << "symbols = " << setprecision(0) << csymbols << endl << endl;
  cout << "samples = " << setprecision(6) << csamples << endl << endl;
  cout << "received symbols = " << setprecision(6) << crsymbols << endl 
    << endl;

  csamples = crc_tx.shape_symbols(csymbols);

  cout << "*** Raised Cosine; complex input ***" << endl << endl;
  cout << "pulse, RC = " << setprecision(6) << crc_tx.get_pulse_shape() 
    << endl << endl;
  cout << "symbols = " << setprecision(0) << csymbols << endl << endl;
  cout << "samples = " << setprecision(6) << csamples << endl << endl;

  return 0;
}
