/*!
 * \file 
 * \brief Convolutional encoder/decoder class test program
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
  cout << "====================================" << endl;
  cout << "    Test of convolutional coders    " << endl;
  cout << "====================================" << endl;

  Array<ivec> spectrum, spectrum_fast, spectrum_punct, spectrum_punct_fast;
  //Array<llvec> spectrum, spectrum_fast, spectrum_punct, spectrum_punct_fast;
  spectrum.set_size(2);
  spectrum_fast.set_size(2);
  spectrum_punct.set_size(2);
  spectrum_punct_fast.set_size(2);

  ivec dist_profile;
  //llvec dist_profile;

  Convolutional_Code code;
  Punctured_Convolutional_Code code_punct;
  BPSK bpsk;

  ivec generator(2);
  generator(0)=0133;
  generator(1)=0171;

  bvec bits=randb(7);
  bvec outbits, temp_outbits, decoded_bits;
  vec tx_symbols, rx_symbols;

  code.set_generator_polynomials(generator, 7);
  code_punct.set_generator_polynomials(generator, 7);

  bmat punct_matrix = "0 1;1 0";
  code_punct.set_puncture_matrix(punct_matrix);

  cout << "code: catastrophic test:" << code.catastrophic() << endl;
  cout << "code: rate:" << code.get_rate() << endl;

  code.calculate_spectrum(spectrum, 10, 10);
  code.fast(spectrum_fast, 10, 10);
  code.distance_profile(dist_profile, 10);

  cout << "spectrum:" << endl;
  cout << "Ad=" << spectrum(0) << endl;
  cout << "Cd=" << spectrum(1) << endl;

  cout << "spectrum, fast:" << endl;
  cout << "Ad=" << spectrum_fast(0) << endl;
  cout << "Cd=" << spectrum_fast(1) << endl;

  cout << "distance profile=:" << dist_profile << endl;

  cout << "code: encode:" << endl;
  code.encode_tail(bits, outbits);

  bpsk.modulate_bits(outbits, tx_symbols);
  rx_symbols = tx_symbols;

  code.decode_tail(tx_symbols, decoded_bits);

  cout << "bits =" << bits << endl;
  cout << "outbits =" << outbits << endl;
  cout << "decoded_bits =" << decoded_bits << endl << endl;

  cout << "====================" << endl;
  cout << "   Punctured code   " << endl;
  cout << "====================" << endl;

  cout << "code_punct: catastrophic test:" << code_punct.catastrophic() << endl;
  cout << "code_punct: rate:" << code_punct.get_rate() << endl;

  cout << "punct_matrix = " << code_punct.get_puncture_matrix() << endl;

  cout << "code_punct: encode:" << endl;
  code_punct.encode_tail(bits, outbits);
  bpsk.modulate_bits(outbits, tx_symbols);
  rx_symbols = tx_symbols;

  code_punct.decode_tail(rx_symbols, decoded_bits);

  cout << "bits =" << bits << endl;
  cout << "outbits =" << outbits << endl;
  cout << "decoded_bits =" << decoded_bits << endl;

  return 0;
}
