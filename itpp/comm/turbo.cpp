/*!
 * \file
 * \brief Implementation of a turbo encoder/decoder class
 * \author Pal Frenger. QLLR support by Erik G. Larsson.
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

#include <itpp/comm/turbo.h>


namespace itpp
{

void Turbo_Codec::set_parameters(ivec gen1, ivec gen2, int constraint_length, const ivec &interleaver_sequence,
                                 int in_iterations, std::string in_metric, double in_logmax_scale_factor,
                                 bool in_adaptive_stop,  LLR_calc_unit in_llrcalc)
{
  //Set the input parameters:
  iterations          = in_iterations;
  interleaver_size    = interleaver_sequence.size();
  Nuncoded            = interleaver_size;
  logmax_scale_factor = in_logmax_scale_factor;
  adaptive_stop       = in_adaptive_stop;

  //Check the decoding metric
  if (in_metric == "LOGMAX") {
    metric = "LOGMAX";
  }
  else if (in_metric == "LOGMAP") {
    metric = "LOGMAP";
  }
  else if (in_metric == "MAP") {
    metric = "MAP";
  }
  else if (in_metric == "TABLE") {
    metric = "TABLE";
  }
  else {
    it_error("Turbo_Codec::set_parameters: The decoder metric must be either MAP, LOGMAP or LOGMAX");
  }

  if (logmax_scale_factor != 1.0) {
    it_assert(metric == "LOGMAX", "Turbo_Codec::set_parameters: logmax_scale_factor can only be used together with LOGMAX decoding");
  }

  //The RSC Encoders:
  rscc1.set_generator_polynomials(gen1, constraint_length);
  rscc2.set_generator_polynomials(gen2, constraint_length);
  n1 = gen1.length() - 1; //Number of parity bits from rscc1
  n2 = gen2.length() - 1; //Number of parity bits from rscc2
  n_tot = 1 + n1 + n2;  //Total number of parity bits and systematic bits

  //Set the number of tail bits:
  m_tail = constraint_length - 1;

  //Calculate the number of coded bits per code-block:
  Ncoded = Nuncoded * n_tot + m_tail * (1 + n1) + m_tail * (1 + n2);

  //Set the interleaver sequence
  bit_interleaver.set_interleaver_depth(interleaver_size);
  float_interleaver.set_interleaver_depth(interleaver_size);
  bit_interleaver.set_interleaver_sequence(interleaver_sequence);
  float_interleaver.set_interleaver_sequence(interleaver_sequence);

  //Default value of the channel reliability scaling factor is 1
  Lc = 1.0;

  // LLR algebra table
  rscc1.set_llrcalc(in_llrcalc);
  rscc2.set_llrcalc(in_llrcalc);

}

void Turbo_Codec::set_interleaver(const ivec &interleaver_sequence)
{
  interleaver_size = interleaver_sequence.size();
  Nuncoded = interleaver_size;

  //Calculate the number of coded bits per code-block:
  Ncoded = Nuncoded * n_tot + m_tail * (1 + n1) + m_tail * (1 + n2);

  //Set the interleaver sequence
  bit_interleaver.set_interleaver_depth(interleaver_size);
  float_interleaver.set_interleaver_depth(interleaver_size);
  bit_interleaver.set_interleaver_sequence(interleaver_sequence);
  float_interleaver.set_interleaver_sequence(interleaver_sequence);
}

void Turbo_Codec::set_metric(std::string in_metric, double in_logmax_scale_factor, LLR_calc_unit in_llrcalc)
{
  logmax_scale_factor = in_logmax_scale_factor;

  //Check the decoding metric
  if (in_metric == "LOGMAX") {
    metric = "LOGMAX";
  }
  else if (in_metric == "LOGMAP") {
    metric = "LOGMAP";
  }
  else if (in_metric == "MAP") {
    metric = "MAP";
  }
  else if (in_metric == "TABLE") {
    metric = "TABLE";
  }
  else {
    it_error("Turbo_Codec::set_metric: The decoder metric must be either MAP, LOGMAP or LOGMAX");
  }

  rscc1.set_llrcalc(in_llrcalc);
  rscc2.set_llrcalc(in_llrcalc);
}

void Turbo_Codec::set_iterations(int in_iterations)
{
  iterations = in_iterations;
}

void Turbo_Codec::set_adaptive_stop(bool in_adaptive_stop)
{
  adaptive_stop = in_adaptive_stop;
}

void Turbo_Codec::set_awgn_channel_parameters(double in_Ec, double in_N0)
{
  Ec = in_Ec;
  N0 = in_N0;
  Lc = 4.0 * std::sqrt(Ec) / N0;
}

void Turbo_Codec::set_scaling_factor(double in_Lc)
{
  Lc = in_Lc;
}


void Turbo_Codec::encode(const bvec &input, bvec &output)
{
  //Local variables:
  int i, k, j, no_blocks;
  int count;
  bvec input_bits, in1, in2, tail1, tail2, out;
  bmat parity1, parity2;

  //Initializations:
  no_blocks = input.length() / Nuncoded;
  output.set_size(no_blocks*Ncoded, false);

  //Set the bit counter to zero:
  count = 0;

  //Encode all code blocks:
  for (i = 0; i < no_blocks; i++) {

    //Encode one block
    input_bits = input.mid(i * Nuncoded, Nuncoded);
    encode_block(input_bits, in1, in2, parity1, parity2);

    //The data part:
    for (k = 0; k < Nuncoded; k++) {
      output(count) = in1(k);
      count++;                                //Systematic bits
      for (j = 0; j < n1; j++) { output(count) = parity1(k, j); count++; } //Parity-1 bits
      for (j = 0; j < n2; j++) { output(count) = parity2(k, j); count++; } //Parity-2 bits
    }

    //The first tail:
    for (k = 0; k < m_tail; k++) {
      output(count) = in1(Nuncoded + k);
      count++;                                //First systematic tail bit
      for (j = 0; j < n1; j++) { output(count) = parity1(Nuncoded + k, j); count++; } //Parity-1 tail bits
    }

    //The second tail:
    for (k = 0; k < m_tail; k++) {
      output(count) = in2(Nuncoded + k);
      count++;                                //Second systematic tail bit
      for (j = 0; j < n2; j++) { output(count) = parity2(Nuncoded + k, j); count++; } //Parity-2 tail bits
    }

  }

}

void Turbo_Codec::decode(const vec &received_signal, bvec &decoded_bits, const bvec &true_bits)
{
  ivec nrof_used_iterations;
  decode(received_signal, decoded_bits, nrof_used_iterations, true_bits);
}

void Turbo_Codec::decode(const vec &received_signal, bvec &decoded_bits, ivec &nrof_used_iterations,
                         const bvec &true_bits)
{

  if ((n1 == 1) && (n2 == 1) && (metric != "MAP")) {
    //This is a speed optimized decoder for R=1/3 (log domain metrics only)
    decode_n3(received_signal, decoded_bits, nrof_used_iterations, true_bits);
  }
  else {

    //Local variables:
    vec rec, rec_syst1, rec_syst2;
    mat rec_parity1, rec_parity2;
    bmat decoded_bits_i;
    int no_blocks, i, j, k, nrof_used_iterations_i;
    int count;
    bool CHECK_TRUE_BITS;

    //Initilaizations:
    no_blocks = received_signal.length() / Ncoded;
    decoded_bits.set_size(no_blocks * Nuncoded, false);
    decoded_bits_i.set_size(iterations, no_blocks * Nuncoded, false);
    rec_syst1.set_size(Nuncoded + m_tail, false);
    rec_syst2.set_size(Nuncoded + m_tail, false);
    rec_syst2.clear();
    rec_parity1.set_size(Nuncoded + m_tail, n1, false);
    rec_parity2.set_size(Nuncoded + m_tail, n2, false);
    nrof_used_iterations.set_size(no_blocks, false);

    //Check the vector true_bits:
    if (true_bits.size() > 1) {
      it_assert(true_bits.size() == (Nuncoded * no_blocks), "Turbo_Codec::decode: Wrong size of input vectors");
      CHECK_TRUE_BITS = true;
    }
    else {
      CHECK_TRUE_BITS = false;
    }

    //Set the bit counter to zero:
    count = 0;

    //Itterate over all received code blocks:
    for (i = 0; i < no_blocks; i++) {

      //The data part:
      for (k = 0; k < Nuncoded; k++) {
        rec_syst1(k) = received_signal(count);
        count++;                               //Systematic bit
        for (j = 0; j < n1; j++) { rec_parity1(k, j) = received_signal(count); count++; }  //Parity-1 bits
        for (j = 0; j < n2; j++) { rec_parity2(k, j) = received_signal(count); count++; }  //Parity-2 bits
      }

      //The first tail:
      for (k = 0; k < m_tail; k++) {
        rec_syst1(Nuncoded + k) = received_signal(count);
        count++;                               //Tail 1 systematic bit
        for (j = 0; j < n1; j++) { rec_parity1(Nuncoded + k, j) = received_signal(count); count++; }  //Tail 1 parity-1 bits
      }

      //The second tail:
      for (k = 0; k < m_tail; k++) {
        rec_syst2(Nuncoded + k) = received_signal(count);
        count++;                              //Tail2 systematic bit
        for (j = 0; j < n2; j++) { rec_parity2(Nuncoded + k, j) = received_signal(count); count++; } //Tali2 parity-2 bits
      }

      //Scale the input data if necessary:
      if (Lc != 1.0) {
        rec_syst1 *= Lc;
        rec_syst2 *= Lc;
        rec_parity1 *= Lc;
        rec_parity2 *= Lc;
      }

      //Decode the block:
      if (CHECK_TRUE_BITS) {
        decode_block(rec_syst1, rec_syst2, rec_parity1, rec_parity2, decoded_bits_i,
                     nrof_used_iterations_i, true_bits.mid(i*Nuncoded, Nuncoded));
        nrof_used_iterations(i) = nrof_used_iterations_i;
      }
      else {
        decode_block(rec_syst1, rec_syst2, rec_parity1, rec_parity2, decoded_bits_i, nrof_used_iterations_i);
        nrof_used_iterations(i) = nrof_used_iterations_i;
      }

      //Put the decoded bits in the output vector:
      decoded_bits.replace_mid(i*Nuncoded, decoded_bits_i.get_row(iterations - 1));

    }

  }

}

void Turbo_Codec::encode_block(const bvec &input, bvec &in1, bvec &in2, bmat &parity1, bmat &parity2)
{
  //Local variables:
  bvec tail1, tail2, interleaved_input;

  //Error check:
  it_assert(input.length() == Nuncoded, "Turbo_Codec::encode_block: Parameter error in Nuncoded.");

  //Initializations:
  tail1.set_size(m_tail, false);
  tail1.clear();
  tail2.set_size(m_tail, false);
  tail2.clear();
  parity1.set_size(Nuncoded + m_tail, n1, false);
  parity1.clear();
  parity2.set_size(Nuncoded + m_tail, n2, false);
  parity2.clear();
  interleaved_input.set_size(Nuncoded, false);
  interleaved_input.clear();

  //The first encoder:
  rscc1.encode_tail(input, tail1, parity1);

  //The interleaver:
  bit_interleaver.interleave(input, interleaved_input);

  //The second encoder:
  rscc2.encode_tail(interleaved_input, tail2, parity2);

  //The input vectors used to the two constituent encoders:
  in1 = concat(input, tail1);
  in2 = concat(interleaved_input, tail2);

}

void Turbo_Codec::decode_block(const vec &rec_syst1, const vec &rec_syst2, const mat &rec_parity1,
                               const mat &rec_parity2, bmat &decoded_bits_i, int &nrof_used_iterations_i,
                               const bvec &true_bits)
{
  //Local variables:
  int i;
  int count, l, k;
  vec extrinsic_input, extrinsic_output, int_rec_syst1, int_rec_syst, tmp;
  vec deint_rec_syst2, rec_syst, sub_rec_syst, Le12, Le21, Le12_int, Le21_int, L, tail1, tail2;
  bool CHECK_TRUE_BITS, CONTINUE;

  //Size initializations:
  decoded_bits_i.set_size(iterations, Nuncoded, false);
  Le12.set_size(Nuncoded + m_tail, false);
  Le21.set_size(Nuncoded + m_tail, false);
  Le21.zeros();

  //Calculate the interleaved and the deinterleaved sequences:
  float_interleaver.interleave(rec_syst1.left(interleaver_size), int_rec_syst1);
  float_interleaver.deinterleave(rec_syst2.left(interleaver_size), deint_rec_syst2);

  //Combine the results from rec_syst1 and rec_syst2 (in case some bits are transmitted several times)
  rec_syst = rec_syst1.left(interleaver_size) + deint_rec_syst2;
  int_rec_syst = rec_syst2.left(interleaver_size) + int_rec_syst1;

  //Get the two tails
  tail1 = rec_syst1.right(m_tail);
  tail2 = rec_syst2.right(m_tail);

  //Form the input vectors (including tails) to the two decoders:
  rec_syst = concat(rec_syst, tail1);
  int_rec_syst = concat(int_rec_syst, tail2);

  // Check the vector true_bits
  if (true_bits.size() > 1) {
    it_assert(true_bits.size() == Nuncoded, "Turbo_Codec::decode_block: Illegal size of input vector true_bits");
    CHECK_TRUE_BITS = true;
  }
  else {
    CHECK_TRUE_BITS = false;
  }

  if (CHECK_TRUE_BITS) {
    it_assert(adaptive_stop == false,
              "Turbo_Codec::decode_block: You can not stop iterations both adaptively and on true bits");
  }

  // Do the iterative decoding:
  nrof_used_iterations_i = iterations;
  for (i = 0; i < iterations; i++) {

    // Decode Code 1
    if (metric == "MAP") {
      rscc1.map_decode(rec_syst, rec_parity1, Le21, Le12, true);
    }
    else if ((metric == "LOGMAX") || (metric == "LOGMAP") || (metric == "TABLE")) {
      rscc1.log_decode(rec_syst, rec_parity1, Le21, Le12, true, metric);
      if (logmax_scale_factor != 1.0) {
        Le12 *= logmax_scale_factor;
      }
    }
    else {
      it_error("Turbo_Codec::decode_block: Illegal metric value");
    }

    // Interleave the extrinsic information:
    float_interleaver.interleave(Le12.left(interleaver_size), tmp);
    Le12_int = concat(tmp, zeros(Le12.size() - interleaver_size));

    // Decode Code 2
    if (metric == "MAP") {
      rscc2.map_decode(int_rec_syst, rec_parity2, Le12_int, Le21_int, true);
    }
    else if ((metric == "LOGMAX") || (metric == "LOGMAP")  || (metric == "TABLE")) {
      rscc2.log_decode(int_rec_syst, rec_parity2, Le12_int, Le21_int, true, metric);
      if (logmax_scale_factor != 1.0) {
        Le21_int *= logmax_scale_factor;
      }
    }
    else {
      it_error("Turbo_Codec::decode_block: Illegal metric value");
    }

    // De-interleave the extrinsic information:
    float_interleaver.deinterleave(Le21_int.left(interleaver_size), tmp);
    Le21 = concat(tmp, zeros(Le21_int.size() - interleaver_size));

    // Take bit decisions
    L = rec_syst + Le21 + Le12;
    count = 0;
    for (l = 0; l < Nuncoded; l++) {
      (L(l) > 0.0) ? (decoded_bits_i(i, count) = bin(0)) : (decoded_bits_i(i, count) = bin(1));
      count++;
    }

    //Check if it is possible to stop iterating early:
    CONTINUE = true;
    if (i < (iterations - 1)) {

      if (CHECK_TRUE_BITS) {
        CONTINUE = false;
        for (k = 0; k < Nuncoded; k++) { if (true_bits(k) != decoded_bits_i(i, k)) { CONTINUE = true; break; } }
      }

      if ((adaptive_stop) && (i > 0)) {
        CONTINUE = false;
        for (k = 0; k < Nuncoded; k++) { if (decoded_bits_i(i - 1, k) != decoded_bits_i(i, k)) { CONTINUE = true; break; } }
      }

    }

    //Check if iterations shall continue:
    if (CONTINUE == false) {
      //Copy the results from current iteration to all following iterations:
      for (k = (i + 1); k < iterations; k++) {
        decoded_bits_i.set_row(k, decoded_bits_i.get_row(i));
        nrof_used_iterations_i = i + 1;
      }
      break;
    }

  }

}

void Turbo_Codec::decode_n3(const vec &received_signal, bvec &decoded_bits, ivec &nrof_used_iterations,
                            const bvec &true_bits)
{
  //Local variables:
  vec rec, rec_syst1, int_rec_syst1, rec_syst2;
  vec rec_parity1, rec_parity2;
  vec extrinsic_input, extrinsic_output, Le12, Le21, Le12_int, Le21_int, L;
  bvec temp_decoded_bits;
  int no_blocks, i, j, k, l, nrof_used_iterations_i;
  int count, count_out;
  bool CHECK_TRUE_BITS, CONTINUE;

  //Initializations:
  no_blocks = received_signal.length() / Ncoded;
  decoded_bits.set_size(no_blocks * Nuncoded, false);
  rec_syst1.set_size(Nuncoded + m_tail, false);
  rec_syst2.set_size(Nuncoded + m_tail, false);
  rec_syst2.clear();
  rec_parity1.set_size(Nuncoded + m_tail, false);
  rec_parity2.set_size(Nuncoded + m_tail, false);
  temp_decoded_bits.set_size(Nuncoded, false);
  decoded_bits_previous_iteration.set_size(Nuncoded, false);
  nrof_used_iterations.set_size(no_blocks, false);

  //Size initializations:
  Le12.set_size(Nuncoded, false);
  Le21.set_size(Nuncoded, false);

  //Set the bit counter to zero:
  count = 0;
  count_out = 0;

  // Check the vector true_bits
  if (true_bits.size() > 1) {
    it_assert(true_bits.size() == Nuncoded*no_blocks, "Turbo_Codec::decode_n3: Illegal size of input vector true_bits");
    CHECK_TRUE_BITS = true;
  }
  else {
    CHECK_TRUE_BITS = false;
  }

  if (CHECK_TRUE_BITS) {
    it_assert(adaptive_stop == false,
              "Turbo_Codec::decode_block: You can not stop iterations both adaptively and on true bits");
  }

  //Iterate over all received code blocks:
  for (i = 0; i < no_blocks; i++) {

    //Reset extrinsic data:
    Le21.zeros();

    //The data part:
    for (k = 0; k < Nuncoded; k++) {
      rec_syst1(k)   = received_signal(count);
      count++; //Systematic bit
      rec_parity1(k) = received_signal(count);
      count++; //Parity-1 bits
      rec_parity2(k) = received_signal(count);
      count++; //Parity-2 bits
    }

    //The first tail:
    for (k = 0; k < m_tail; k++) {
      rec_syst1(Nuncoded + k)   = received_signal(count);
      count++; //Tail 1 systematic bit
      rec_parity1(Nuncoded + k) = received_signal(count);
      count++; //Tail 1 parity-1 bits
    }

    //The second tail:
    for (k = 0; k < m_tail; k++) {
      rec_syst2(Nuncoded + k)   = received_signal(count);
      count++; //Tail2 systematic bit
      rec_parity2(Nuncoded + k) = received_signal(count);
      count++; //Tali2 parity-2 bits
    }

    float_interleaver.interleave(rec_syst1.left(Nuncoded), int_rec_syst1);
    rec_syst2.replace_mid(0, int_rec_syst1);

    //Scale the input data if necessary:
    if (Lc != 1.0) {
      rec_syst1   *= Lc;
      rec_syst2   *= Lc;
      rec_parity1 *= Lc;
      rec_parity2 *= Lc;
    }

    //Decode the block:
    CONTINUE = true;
    nrof_used_iterations_i = iterations;
    for (j = 0; j < iterations; j++) {

      rscc1.log_decode_n2(rec_syst1, rec_parity1, Le21, Le12, true, metric);
      if (logmax_scale_factor != 1.0) { Le12 *= logmax_scale_factor; }
      float_interleaver.interleave(Le12, Le12_int);

      rscc2.log_decode_n2(rec_syst2, rec_parity2, Le12_int, Le21_int, true, metric);
      if (logmax_scale_factor != 1.0) { Le21_int *= logmax_scale_factor; }
      float_interleaver.deinterleave(Le21_int, Le21);

      if (adaptive_stop) {
        L = rec_syst1.left(Nuncoded) + Le21.left(Nuncoded) + Le12.left(Nuncoded);
        for (l = 0; l < Nuncoded; l++) {(L(l) > 0.0) ? (temp_decoded_bits(l) = bin(0)) : (temp_decoded_bits(l) = bin(1)); }
        if (j == 0) { decoded_bits_previous_iteration = temp_decoded_bits; }
        else {
          if (temp_decoded_bits == decoded_bits_previous_iteration) {
            CONTINUE = false;
          }
          else if (j < (iterations - 1)) {
            decoded_bits_previous_iteration = temp_decoded_bits;
          }
        }
      }

      if (CHECK_TRUE_BITS) {
        L = rec_syst1.left(Nuncoded) + Le21.left(Nuncoded) + Le12.left(Nuncoded);
        for (l = 0; l < Nuncoded; l++) {(L(l) > 0.0) ? (temp_decoded_bits(l) = bin(0)) : (temp_decoded_bits(l) = bin(1)); }
        if (temp_decoded_bits == true_bits.mid(i*Nuncoded, Nuncoded)) {
          CONTINUE = false;
        }
      }

      if (CONTINUE == false) { nrof_used_iterations_i = j + 1; break; }

    }

    //Take final bit decisions
    L = rec_syst1.left(Nuncoded) + Le21.left(Nuncoded) + Le12.left(Nuncoded);
    for (l = 0; l < Nuncoded; l++) {
      (L(l) > 0.0) ? (decoded_bits(count_out) = bin(0)) : (decoded_bits(count_out) = bin(1));
      count_out++;
    }

    nrof_used_iterations(i) = nrof_used_iterations_i;

  }

}

ivec wcdma_turbo_interleaver_sequence(int interleaver_size)
{
  const int MAX_INTERLEAVER_SIZE = 5114;
  const int MIN_INTERLEAVER_SIZE = 40;
  int K;  //Interleaver size
  int R;  //Number of rows of rectangular matrix
  int C;  //Number of columns of rectangular matrix
  int p;  //Prime number
  int v;  //Primitive root
  ivec s; //Base sequence for intra-row permutation
  ivec q; //Minimum prime integers
  ivec r; //Permuted prime integers
  ivec T; //Inter-row permutation pattern
  imat U; //Intra-row permutation patter
  ivec I; //The interleaver sequence
  ivec primes, roots, Pat1, Pat2, Pat3, Pat4, Isort;
  int i, j, qj, temp, row, col, index, count;

  if (interleaver_size > MAX_INTERLEAVER_SIZE) {

    I = sort_index(randu(interleaver_size));
    return I;

  }
  else {

    p = 0;
    v = 0;

    //Check the range of the interleaver size:
    it_assert(interleaver_size <= MAX_INTERLEAVER_SIZE, "wcdma_turbo_interleaver_sequence: The interleaver size is to large");
    it_assert(interleaver_size >= MIN_INTERLEAVER_SIZE, "wcdma_turbo_interleaver_sequence: The interleaver size is to small");

    K = interleaver_size;

    //Definitions of primes and associated primitive roots:
    primes = "2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 233 239 241 251 257";
    roots = "0 0 0 3 2 2 3 2 5 2 3 2 6 3 5 2 2 2 2 7 5 3 2 3 5 2 5 2 6 3 3 2 3 2 2 6 5 2 5 2 2 2 19 5 2 3 2 3 2 6 3 7 7 6 3";

    //Determine R
    if ((K >= 40) && (K <= 159)) {
      R = 5;
    }
    else if (((K >= 160) && (K <= 200)) || ((K >= 481) && (K <= 530))) {
      R = 10;
    }
    else {
      R = 20;
    }

    //Determine C
    if ((K >= 481) && (K <= 530)) {
      p = 53;
      v = 2;
      C = p;
    }
    else {
      //Find minimum prime p such that (p+1) - K/R >= 0 ...
      for (i = 0; i < primes.length(); i++) {
        if ((double(primes(i) + 1) - double(K) / double(R)) >= 0.0) {
          p = primes(i);
          v = roots(i);
          break;
        }
      }
      //... and etermine C such that
      if ((double(p) - double(K) / double(R)) >= 0.0) {
        if ((double(p) - 1.0 - double(K) / double(R)) >= 0.0) {
          C = p - 1;
        }
        else {
          C = p;
        }
      }
      else {
        C = p + 1;
      }
    }

    //Construct the base sequencs s for intra-row permutaions
    s.set_size(p - 1, false);
    s.clear();
    s(0) = 1;
    for (i = 1; i <= (p - 2); i++) {
      s(i) = mod(v * s(i - 1), p);
    }

    //Let q(0) = 1 be the first prime integer in {q(j)}, and select the consecutive
    //minimum prime integers {q(j)}, j = 1, 2, ..., (R-1) such that gcd( q(j), p-1) == 1, q(j) > 6, and q(j) > q(j-1)
    q.set_size(R, false);
    q.clear();
    q(0) = 1;
    for (j = 1; j <= (R - 1); j++) {
      for (i = 0; i < primes.length(); i++) {
        qj = primes(i);
        if ((qj > 6) && (qj > q(j - 1))) {
          if (gcd(qj, p - 1) == 1) {
            q(j) = qj;
            break;
          }
        }
      }
    }

    //Definitions of Pat1, Pat2, Pat3, and Pat4:
    Pat1 = "19 9 14 4 0 2 5 7 12 18 10 8 13 17 3 1 16 6 15 11";
    Pat2 = "19 9 14 4 0 2 5 7 12 18 16 13 17 15 3 1 6 11 8 10";
    Pat3 = "9 8 7 6 5 4 3 2 1 0";
    Pat4 = "4 3 2 1 0";

    //T(j) is the inter-row permutation patters defined as one of the following four
    //kinds of patterns: Pat1, Pat2, Pat3, and Pat4 depending on the number of input bits K
    if (K >= 3211) {
      T = Pat1;
    }
    else if (K >= 3161) {
      T = Pat2;
    }
    else if (K >= 2481) {
      T = Pat1;
    }
    else if (K >= 2281) {
      T = Pat2;
    }
    else if (K >= 531) {
      T = Pat1;
    }
    else if (K >= 481) {
      T = Pat3;
    }
    else if (K >= 201) {
      T = Pat1;
    }
    else if (K >= 160) {
      T = Pat3;
    }
    else {
      T = Pat4;
    }

    //Permute {q(j)} to make {r(j)} such that r(T(j)) = q(j), j = 0, 1, ..., (R-1),
    //where T(j) indicates the original row position of the j-th permuted row
    r.set_size(R, false);
    r.clear();
    for (j = 0; j <= (R - 1); j++) {
      r(T(j)) = q(j);
    }

    //U(j,i) is the input bit position of i-th output after the permutation of j-th row
    //Perform the j-th (j=0, 1, 2, ..., (R-1)) intra-row permutation as
    U.set_size(R, C, false);
    U.clear();
    if (C == p) {
      for (j = 0; j <= (R - 1); j++) {
        for (i = 0; i <= (p - 2); i++) {
          U(j, i) = s(mod(i * r(j), p - 1));
        }
        U(j, p - 1) = 0;
      }
    }
    else if (C == (p + 1)) {
      for (j = 0; j <= (R - 1); j++) {
        for (i = 0; i <= (p - 2); i++) {
          U(j, i) = s(mod(i * r(j), p - 1));
        }
        U(j, p - 1) = 0;
        U(j, p) = p;
      }
      if (K == (C*R)) {
        temp = U(R - 1, p);
        U(R - 1, p) = U(R - 1, 0);
        U(R - 1, 0) = temp;
      }
    }
    else if (C == (p - 1)) {
      for (j = 0; j <= (R - 1); j++) {
        for (i = 0; i <= (p - 2); i++) {
          U(j, i) = s(mod(i * r(j), p - 1)) - 1;
        }
      }
    }

    //Calculate the interleaver sequence:
    I.set_size(K, false);
    I.clear();
    count = 0;
    for (i = 0; i < C; i++) {
      for (j = 0; j < R; j++) {
        row = T(j);
        col = U(row, i);
        index = row * C + col;
        if (index < K) {
          I(count) = index;
          count++;
        }
      }
    }

    return I;
  }
}

} // namespace itpp
