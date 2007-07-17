/*!
 * \file
 * \brief 1D and 2D modulators test program
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;


int main()
{
  RNG_reset(12345);

  cout << "===========================================================" << endl;
  cout << "                    Test of Modulators                     " << endl;
  cout << "===========================================================" << endl;

  const int no_symbols = 5;
  const double N0 = 0.1;

  {
    cout << endl << "Modulator_1D (configured as BPSK)" << endl;
    Modulator_1D  mod("1.0 -1.0", "0 1");
    int bps = round_i(mod.bits_per_symbol());

    bvec tx_bits = randb(no_symbols * bps);
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);
    vec noise = sqrt(N0) * randn(no_symbols);

    vec tx_symbols = mod.modulate_bits(tx_bits);
    vec rx_symbols = tx_symbols + noise;
    bvec decbits = mod.demodulate_bits(rx_symbols);

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  decbits         = " << decbits << endl;

    tx_symbols = mod.modulate(tx_sym_numbers);
    rx_symbols = tx_symbols + noise;
    ivec dec_sym_numbers = mod.demodulate(rx_symbols);

    cout << "* modulating symbol numbers:" << endl;
    cout << "  tx_sym_numbers  = " << tx_sym_numbers << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  dec_sym_numbers = " << dec_sym_numbers << endl;

    cout << endl << "BPSK (real signal)" << endl;
    BPSK bpsk;

    bpsk.modulate_bits(tx_bits, tx_symbols);
    rx_symbols = tx_symbols + noise;
    bpsk.demodulate_bits(rx_symbols, decbits);

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  decbits         = " << decbits << endl;

    cout << endl << "BPSK (complex signal)" << endl;
    BPSK_c bpsk_c;

    cvec tx_csymbols = bpsk_c.modulate_bits(tx_bits);
    cvec rx_csymbols = tx_csymbols + to_cvec(noise, -noise);
    decbits = bpsk_c.demodulate_bits(rx_csymbols);
    vec softbits_approx = bpsk_c.demodulate_soft_bits(rx_csymbols, N0, APPROX);
    vec softbits = bpsk_c.demodulate_soft_bits(rx_csymbols, N0, LOGMAP);

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_csymbols     = " << tx_csymbols << endl;
    cout << "  rx_csymbols     = " << rx_csymbols << endl;
    cout << "  decbits         = " << decbits << endl;
    cout << "  softbits        = " << softbits << endl;
    cout << "  softbits_approx = " << softbits_approx << endl << endl;
  }

  cout << "===========================================================" << endl;

  {
    cout << endl << "Modulator_1D (configured as 4-PAM)" << endl;
    Modulator_1D mod("-3.0 -1.0 1.0 3.0", "0 1 3 2");
    int bps = round_i(mod.bits_per_symbol());

    bvec tx_bits = randb(no_symbols * bps);
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);
    vec noise = sqrt(N0) * randn(no_symbols);

    vec tx_symbols = mod.modulate_bits(tx_bits);
    vec rx_symbols = tx_symbols + noise;
    bvec decbits = mod.demodulate_bits(rx_symbols);

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  decbits         = " << decbits << endl;

    tx_symbols = mod.modulate(tx_sym_numbers);
    rx_symbols = tx_symbols + noise;
    ivec dec_sym_numbers = mod.demodulate(rx_symbols);

    cout << "* modulating symbol numbers:" << endl;
    cout << "  tx_sym_numbers  = " << tx_sym_numbers << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  dec_sym_numbers = " << dec_sym_numbers << endl;

    cout << endl << "4-PAM (real signal)" << endl;
    PAM pam(4);

    pam.modulate_bits(tx_bits, tx_symbols);
    rx_symbols = tx_symbols + noise;
    pam.demodulate_bits(rx_symbols, decbits);

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  decbits         = " << decbits << endl;

    cout << endl << "4-PAM (complex signal)" << endl;
    PAM_c pam_c(4);

    cvec tx_csymbols = pam_c.modulate_bits(tx_bits);
    cvec rx_csymbols = tx_csymbols + to_cvec(noise, -noise);
    decbits = pam_c.demodulate_bits(rx_csymbols);
    vec softbits_approx = pam_c.demodulate_soft_bits(rx_csymbols, N0, APPROX);
    vec softbits = pam_c.demodulate_soft_bits(rx_csymbols, N0, LOGMAP);

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_csymbols     = " << tx_csymbols << endl;
    cout << "  rx_csymbols     = " << rx_csymbols << endl;
    cout << "  decbits         = " << decbits << endl;
    cout << "  softbits        = " << softbits << endl;
    cout << "  softbits_approx = " << softbits_approx << endl << endl;
  }

  cout << "===========================================================" << endl;

  {
    cout << endl << "Modulator_2D (configured as QPSK)" << endl;
    Modulator_2D mod("1+0i 0+1i -1+0i 0-1i", "0 1 3 2");
    int bps = round_i(mod.bits_per_symbol());

    bvec tx_bits = randb(no_symbols * bps);
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);
    cvec noise = sqrt(N0) * randn_c(no_symbols);

    cvec tx_symbols = mod.modulate(tx_sym_numbers);
    cvec rx_symbols = tx_symbols + noise;
    ivec dec_sym_numbers = mod.demodulate(rx_symbols);

    cout << "* modulating symbol numbers:" << endl;
    cout << "  tx_sym_numbers  = " << tx_sym_numbers << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  dec_sym_numbers = " << dec_sym_numbers << endl;

    tx_symbols = mod.modulate_bits(tx_bits);
    rx_symbols = tx_symbols + noise;
    bvec decbits = mod.demodulate_bits(rx_symbols);
    vec softbits_approx = mod.demodulate_soft_bits(rx_symbols, N0, APPROX);
    vec softbits = mod.demodulate_soft_bits(rx_symbols, N0, LOGMAP);

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  decbits         = " << decbits << endl;
    cout << "  softbits        = " << softbits << endl;
    cout << "  softbits_approx = " << softbits_approx << endl;

    cout << endl << "QPSK" << endl;
    QPSK qpsk;

    tx_symbols = qpsk.modulate(tx_sym_numbers);
    rx_symbols = tx_symbols + noise;
    dec_sym_numbers = qpsk.demodulate(rx_symbols);

    cout << "* modulating symbol numbers:" << endl;
    cout << "  tx_sym_numbers  = " << tx_sym_numbers << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  dec_sym_numbers = " << dec_sym_numbers << endl;

    tx_symbols = qpsk.modulate_bits(tx_bits);
    rx_symbols = tx_symbols + noise;
    decbits = qpsk.demodulate_bits(rx_symbols);
    softbits_approx = qpsk.demodulate_soft_bits(rx_symbols, N0, APPROX);
    softbits = qpsk.demodulate_soft_bits(rx_symbols, N0, LOGMAP);

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  decbits         = " << decbits << endl;
    cout << "  softbits        = " << softbits << endl;
    cout << "  softbits_approx = " << softbits_approx << endl << endl;
  }

  cout << "===========================================================" << endl;

  {
    cout << endl << "8-PSK" << endl;
    PSK psk(8);
    int bps = round_i(psk.bits_per_symbol());

    bvec tx_bits = randb(no_symbols * bps);
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);
    cvec noise = sqrt(N0) * randn_c(no_symbols);

    cvec tx_symbols = psk.modulate(tx_sym_numbers);
    cvec rx_symbols = tx_symbols + noise;
    ivec dec_sym_numbers = psk.demodulate(rx_symbols);

    cout << "* modulating symbol numbers:" << endl;
    cout << "  tx_sym_numbers  = " << tx_sym_numbers << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  dec_sym_numbers = " << dec_sym_numbers << endl;

    tx_symbols = psk.modulate_bits(tx_bits);
    rx_symbols = tx_symbols + noise;
    bvec decbits = psk.demodulate_bits(rx_symbols);
    vec softbits_approx = psk.demodulate_soft_bits(rx_symbols, N0, APPROX);
    vec softbits = psk.demodulate_soft_bits(rx_symbols, N0, LOGMAP);

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  decbits         = " << decbits << endl;
    cout << "  softbits        = " << softbits << endl;
    cout << "  softbits_approx = " << softbits_approx << endl << endl;
  }

  cout << "===========================================================" << endl;

  {
    cout << endl << "16-QAM" << endl;
    QAM qam(16);
    int bps = round_i(qam.bits_per_symbol());

    bvec tx_bits = randb(no_symbols * bps);
    ivec tx_sym_numbers = randi(no_symbols, 0, pow2i(bps) - 1);
    cvec noise = sqrt(N0) * randn_c(no_symbols);

    cvec tx_symbols = qam.modulate(tx_sym_numbers);
    cvec rx_symbols = tx_symbols + noise;
    ivec dec_sym_numbers = qam.demodulate(rx_symbols);

    cout << "* modulating symbol numbers:" << endl;
    cout << "  tx_sym_numbers  = " << tx_sym_numbers << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  dec_sym_numbers = " << dec_sym_numbers << endl;

    tx_symbols = qam.modulate_bits(tx_bits);
    rx_symbols = tx_symbols + noise;
    bvec decbits = qam.demodulate_bits(rx_symbols);
    vec softbits_approx = qam.demodulate_soft_bits(rx_symbols, N0, APPROX);
    vec softbits = qam.demodulate_soft_bits(rx_symbols, N0, LOGMAP);

    cout << "* modulating bits:" << endl;
    cout << "  tx_bits         = " << tx_bits << endl;
    cout << "  tx_symbols      = " << tx_symbols << endl;
    cout << "  rx_symbols      = " << rx_symbols << endl;
    cout << "  decbits         = " << decbits << endl;
    cout << "  softbits        = " << softbits << endl;
    cout << "  softbits_approx = " << softbits_approx << endl << endl;
  }
}
