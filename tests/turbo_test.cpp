/*!
 * \file 
 * \brief Turbo encoder/decoder class test program
 * \author Pal Frenger
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

#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;


int main()
{
    Turbo_Codec turbo;
    ivec gen(2);
    gen(0) = 07; gen(1) = 05;
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
    vec ber(EbN0db.length()); ber.clear();
    vec avg_nrof_iterations(EbN0db.length()); avg_nrof_iterations.clear();
    ivec nrof_used_iterations;
    int i;

    vec rec_systematic, rec_parity, rec_parity1, rec_parity2, extrinsic_input, 
      extrinsic_output, received, symbols;
    bvec input, parity1, parity2, parity_punct, tail1, output1, output2, 
      transmitted, uncoded_bits;
    bvec coded_bits, decoded_bits;

    Normal_RNG noise_src;
    RNG_reset(12345);

    BPSK bpsk;
    BERC berc;

    cout << "=============================================" << endl;
    cout << "           Starting Simulation               " << endl;
    cout << " Bit error rate as a function of itterations " << endl;
    cout << "=============================================" << endl;
    cout << "  Block length = " << block_length << endl;
    cout << "  Generator polynomials = " << std::oct << gen << std::dec << endl;
    cout << "  Max number of Iterations = " << iterations << endl;
    cout << "  Adaptive stop = " << adaptive_stop << endl;
    cout << "  Eb/N0 = " << EbN0db << " [dB]" << endl;
    cout << "  Turbo encoder rate 1/3 (plus tail bits)" << endl;
    cout << "=============================================" << endl;

    vec err(EbN0db.length());
    vec cor(EbN0db.length());

    for (i = 0; i < EbN0db.length(); i++) {
	cout << "Now simulating EbN0db = " << EbN0db(i) << endl;

	noise_src.setup(0.0, sigma2(i));
	turbo.set_awgn_channel_parameters(Ec, N0(i));
	input = randb(block_length * num_blocks);

	turbo.encode(input,transmitted);
	bpsk.modulate_bits(transmitted, symbols);
	received = symbols + noise_src(transmitted.length());
	turbo.decode(received, decoded_bits, nrof_used_iterations);

	berc.clear();
	berc.count(input,decoded_bits);
	err(i) = berc.get_errors();
	cor(i) = berc.get_corrects();
	ber(i) = berc.get_errorrate();
	avg_nrof_iterations(i) = static_cast<double>(sum(nrof_used_iterations)) 
	  / length(nrof_used_iterations);
    }

    cout << "Average numer of iterations used = " << avg_nrof_iterations 
	 << endl;
    cout << "err = " << err << endl;
    cout << "cor = " << cor << endl;
    cout << "ber = " << ber << endl;

    return 0;
}
