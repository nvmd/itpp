/*!
 * \file
 * \brief Turbo encoder/decoder class test program
 * \author Pal Frenger and Erik G. Larsson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;


int main()
{
  Turbo_Codec turbo;
  ivec gen(2);
  gen(0) = 07;
  gen(1) = 05;
  int constraint_length = 3;
  int block_length = 400;
  ivec interleaver_sequence = wcdma_turbo_interleaver_sequence(block_length);
  int iterations = 8;
  string metric = "LOGMAX";
  double logmax_scale_factor = 0.7;
  bool adaptive_stop = true;
  turbo.set_parameters(gen, gen, constraint_length, interleaver_sequence,
                       iterations, metric, logmax_scale_factor,
                       adaptive_stop);
  int num_blocks = 50;

  vec EbN0db = "0.0 0.5 1.0 1.5 2.0";
  double A = 1.0;
  double Ts = 1.0;
  double Ec = A * A * Ts;
  double r = block_length / (3.0 * block_length + 8.0); // 8 tailbits
  double Eb = Ec / r;
  vec EbN0 = pow(10.0, 0.1 * EbN0db);
  vec N0 = Eb * pow(EbN0, -1.0);
  vec sigma2 = N0 / 2;
  ivec nrof_used_iterations;

  vec symbols, received;
  bvec input, coded_bits, decoded_bits, transmitted;

  Normal_RNG noise_src;
  RNG_reset(12345);

  BPSK bpsk;
  BERC berc;

  cout << "=============================================" << endl;
  cout << "           Starting Simulation               " << endl;
  cout << " Bit error rate as a function of Eb/N0       " << endl;
  cout << "=============================================" << endl;
  cout << "  Block length = " << block_length << endl;
  cout << "  Generator polynomials = " << std::oct << gen << std::dec << endl;
  cout << "  Max number of Iterations = " << iterations << endl;
  cout << "  Adaptive stop = " << adaptive_stop << endl;
  cout << "  Eb/N0 = " << EbN0db << " [dB]" << endl;
  cout << "  Turbo encoder rate 1/3 (plus tail bits)" << endl;
  cout << "=============================================" << endl;

  mat err = zeros(4, EbN0db.length());
  mat cor = zeros(4, EbN0db.length());
  mat ber = zeros(4, EbN0db.length());
  mat avg_nrof_iterations = zeros(4, EbN0db.length());
  LLR_calc_unit lowresllrcalc(10, 7, 9);  // table with low resolution
  Array<Real_Timer> timer(4);
  for (int i = 0; i < 4; i++) {
    timer(i).reset();
  }

  for (int i = 0; i < EbN0db.length(); i++) {
    cout << "Now simulating EbN0db = " << EbN0db(i) << endl;

    noise_src.setup(0.0, sigma2(i));
    turbo.set_awgn_channel_parameters(Ec, N0(i));
    input = randb(block_length * num_blocks);

    turbo.encode(input, transmitted);
    bpsk.modulate_bits(transmitted, symbols);
    received = symbols + noise_src(transmitted.length());

    // -- logmax decoding --
    turbo.set_metric("LOGMAX", 1.0);
    timer(0).start();
    turbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(0).stop();
    berc.clear();
    berc.count(input, decoded_bits);
    err(0, i) = berc.get_errors();
    cor(0, i) = berc.get_corrects();
    ber(0, i) = berc.get_errorrate();
    avg_nrof_iterations(0, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

    // -- logmap decoding --
    turbo.set_metric("LOGMAP", 1.0);
    timer(1).start();
    turbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(1).stop();
    berc.clear();
    berc.count(input, decoded_bits);
    err(1, i) = berc.get_errors();
    cor(1, i) = berc.get_corrects();
    ber(1, i) = berc.get_errorrate();
    avg_nrof_iterations(1, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

    // -- QLLR decoding, default resolution --
    turbo.set_metric("TABLE", 1.0);
    timer(2).start();
    turbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(2).stop();
    berc.clear();
    berc.count(input, decoded_bits);
    err(2, i) = berc.get_errors();
    cor(2, i) = berc.get_corrects();
    ber(2, i) = berc.get_errorrate();
    avg_nrof_iterations(2, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

    // -- QLLR decoding, low resolution --
    turbo.set_metric("TABLE", 1.0, lowresllrcalc);
    timer(3).start();
    turbo.decode(received, decoded_bits, nrof_used_iterations);
    timer(3).stop();
    berc.clear();
    berc.count(input, decoded_bits);
    err(3, i) = berc.get_errors();
    cor(3, i) = berc.get_corrects();
    ber(3, i) = berc.get_errorrate();
    avg_nrof_iterations(3, i) = static_cast<double>(sum(nrof_used_iterations)) / length(nrof_used_iterations);

  }

  cout << "Results: (1st row: logmax, 2nd row: logmap, 3rd row: qllr, default resolution, 4th row: qllr, low resolution" << endl;
  cout << "Bit error rate: " << endl;
  cout << "ber = " << ber << endl;
  cout << "Average numer of iterations used: " << endl;
  cout << avg_nrof_iterations  << endl;
  cout << "Number of bit errors counted: " << endl;
  cout << "err = " << err << endl;
  cout << "Number of correct bits counted: " << endl;
  cout << "cor = " << cor << endl;

  /*
  // The test program cannot print this out, but on my system
  // the QLLR based decoder is about 8 times faster than logmap. -EGL
  cout << "Timers: ";
  for (int i=0; i<4; i++) { cout << timer(i).get_time() << "  "; }
  cout << endl;
  */

  return 0;
}
