/*!
 * \file 
 * \brief Modulator/demodulator classes test program
 * \author Tony Ottosson and Adam Piatyszek
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * -------------------------------------------------------------------------
 */

#include <itpp/itbase.h>
#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;


int main()
{
  RNG_reset(12345);

  cout << "===========================================================" << endl;
  cout << "                    Test of Modulators                     " << endl;
  cout << "===========================================================" << endl;
   
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






