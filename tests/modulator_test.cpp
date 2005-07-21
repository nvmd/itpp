#include <itpp/itcomm.h>

using std::cout;
using std::endl;
using namespace itpp;

int main()
{
  RNG_reset(12345);

    cout << "=============================================================" << endl;
    cout << "                    Test of Modulators                       " << endl;
    cout << "=============================================================" << endl;
    
    int no_bits = 5;

  {
    cout << endl << "BPSK" << endl;
    BPSK mod;
    
    bvec tx_bits = randb(no_bits);
    cvec tx_symbols = mod.modulate_bits(tx_bits);
    double N0 = 0.1;
    cvec noise = sqrt(N0/2.0) * randn_c(tx_symbols.size());
    cvec rx_symbols = tx_symbols + noise;
    
    bvec decbits;
    vec softbits;
    mod.demodulate_bits(rx_symbols,decbits);
    mod.demodulate_soft_bits(rx_symbols,N0,softbits);
    
    cout << "tx_bits         = " << tx_bits << endl;
    cout << "tx_symbols      = " << tx_symbols << endl;
    cout << "rx_symbols      = " << rx_symbols << endl;
    cout << "decbits         = " << decbits << endl;
    cout << "softbits        = " << softbits << endl;
  }
  {
    cout << endl << "QPSK" << endl;
    QPSK mod;
    
    bvec tx_bits = randb(no_bits*2);
    cvec tx_symbols = mod.modulate_bits(tx_bits);
    double N0 = 0.1;
    cvec noise = sqrt(N0/2.0) * randn_c(tx_symbols.size());
    cvec rx_symbols = tx_symbols + noise;
    
    bvec decbits;
    vec softbits;
    mod.demodulate_bits(rx_symbols,decbits);
    mod.demodulate_soft_bits(rx_symbols,N0,softbits);
    
    cout << "tx_bits         = " << tx_bits << endl;
    cout << "tx_symbols      = " << tx_symbols << endl;
    cout << "rx_symbols      = " << rx_symbols << endl;
    cout << "decbits         = " << decbits << endl;
    cout << "softbits        = " << softbits << endl;
  }
  {
    cout << endl << "8-PSK" << endl;
    PSK mod(8);
    
    bvec tx_bits = randb(no_bits*3);
    cvec tx_symbols = mod.modulate_bits(tx_bits);
    double N0 = 0.1;
    cvec noise = sqrt(N0/2.0) * randn_c(tx_symbols.size());
    cvec rx_symbols = tx_symbols + noise;
    
    bvec decbits;
    vec softbits;
    mod.demodulate_bits(rx_symbols,decbits);
    mod.demodulate_soft_bits(rx_symbols,N0,softbits);
    
    cout << "tx_bits         = " << tx_bits << endl;
    cout << "tx_symbols      = " << tx_symbols << endl;
    cout << "rx_symbols      = " << rx_symbols << endl;
    cout << "decbits         = " << decbits << endl;
    cout << "softbits        = " << softbits << endl;
  }
  {
    cout << endl << "8-PAM" << endl;
    PAM mod(8);
    
    bvec tx_bits = randb(no_bits*3);
    cvec tx_symbols = mod.modulate_bits(tx_bits);
    double N0 = 0.1;
    cvec noise = sqrt(N0/2.0) * randn_c(tx_symbols.size());
    cvec rx_symbols = tx_symbols + noise;
    
    bvec decbits;
    vec softbits;
    mod.demodulate_bits(rx_symbols,decbits);
    mod.demodulate_soft_bits(rx_symbols,N0,softbits);
    
    cout << "tx_bits         = " << tx_bits << endl;
    cout << "tx_symbols      = " << tx_symbols << endl;
    cout << "rx_symbols      = " << rx_symbols << endl;
    cout << "decbits         = " << decbits << endl;
    cout << "softbits        = " << softbits << endl;
  }
  {
    cout << endl << "16-QAM" << endl;
    QAM mod(16);
    
    bvec tx_bits = randb(no_bits*4);
    cvec tx_symbols = mod.modulate_bits(tx_bits);
    double N0 = 0.1;
    cvec noise = sqrt(N0/2.0) * randn_c(tx_symbols.size());
    cvec rx_symbols = tx_symbols + noise;
    
    bvec decbits;
    vec softbits;
    mod.demodulate_bits(rx_symbols,decbits);
    mod.demodulate_soft_bits(rx_symbols,N0,softbits);
    
    cout << "tx_bits         = " << tx_bits << endl;
    cout << "tx_symbols      = " << tx_symbols << endl;
    cout << "rx_symbols      = " << rx_symbols << endl;
    cout << "decbits         = " << decbits << endl;
    cout << "softbits        = " << softbits << endl;
  }
}






