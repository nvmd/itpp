#include "itpp/itcomm.h"

using std::cout;
using std::endl;
using std::complex;
using namespace itpp;

int main()
{
  cout << "======================================" << endl;
  cout << "    Test of pulse shaping routines" << endl;
  cout << "======================================" << endl;

  BPSK bpsk;
  vec symbols, samples, rsymbols;
  Root_Raised_Cosine<double> rrc_tx(0.5), rrc_rx(0.5);
  Raised_Cosine<double> rc_tx(0.5);

  bpsk.modulate_bits(randb(20), symbols);
  samples = rrc_tx.shape_symbols(symbols);
  rsymbols = rrc_rx.shape_samples(samples);

  cout << "============== Root Raised Cosine ===================" << endl;
  cout << "pulse, RRC=" << rrc_tx.get_pulse_shape() << endl;
  cout << "symbols=" << symbols << endl;
  cout << "samples=" << samples << endl;
  cout << "received symbols=" << rsymbols << endl;

  samples = rc_tx.shape_symbols(symbols);

  cout << "============== Raised Cosine ===================" << endl;
  cout << "pulse, RC=" << rc_tx.get_pulse_shape() << endl;
  cout << "symbols=" << symbols << endl;
  cout << "samples=" << samples << endl;


  QPSK qpsk;
  cvec csymbols, csamples, crsymbols;
  Root_Raised_Cosine<complex<double> > crrc_tx(0.5), crrc_rx(0.5);
  Raised_Cosine<complex<double> > crc_tx(0.5);


  qpsk.modulate_bits(randb(40), csymbols);
  csamples = crrc_tx.shape_symbols(csymbols);
  crsymbols = crrc_rx.shape_samples(csamples);

  cout << "============== Root Raised Cosine ===================" << endl;
  cout << "pulse, RRC=" << crrc_tx.get_pulse_shape() << endl;
  cout << "symbols=" << csymbols << endl;
  cout << "samples=" << csamples << endl;
  cout << "received symbols=" << crsymbols << endl;


  csamples = crc_tx.shape_symbols(csymbols);

  cout << "============== Raised Cosine ===================" << endl;
  cout << "pulse, RC=" << crc_tx.get_pulse_shape() << endl;
  cout << "symbols=" << csymbols << endl;
  cout << "samples=" << csamples << endl;

  return 0;
}
