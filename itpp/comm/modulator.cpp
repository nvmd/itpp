/*!
 * \file 
 * \brief Implementation of modulator classes
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
 * Copyright (C) 1995-2006  (see AUTHORS file for a list of contributors)
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

#include <itpp/comm/modulator.h>
#include <itpp/comm/commfunc.h>
#include <itpp/base/converters.h>
#include <itpp/base/logexpfunc.h>
#include <itpp/base/elmatfunc.h>

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

namespace itpp {

  //------------- class: Modulator_1d ----------------

  Modulator_1d::Modulator_1d(const vec &symbols, const ivec &bitmap)
  {
    set(symbols, bitmap);
  }

  void Modulator_1d::set(const vec &in_symbols, const ivec &in_bitmap)
  {
    it_assert(in_symbols.size() == in_bitmap.size(), "Modulator_1d::set(): Number of symbols and bitmap does not match");
    M = in_bitmap.size();
    k = needed_bits(M - 1);
    symbols = in_symbols;
    bitmap = in_bitmap; 
  }

  vec Modulator_1d::modulate(const ivec &symbolnumbers) const
  {
    vec temp(symbolnumbers.size());
    for (int i = 0; i < symbolnumbers.size(); i++)
      temp(i) = symbols(symbolnumbers(i));
    return temp;
  }

  vec Modulator_1d::modulate_bits(const bvec &bits) const
  {
    int no_symbols = bits.length() / k;
    int pos, symb;
    vec output = zeros(no_symbols);

    for (int i = 0; i < no_symbols; i++) {
      pos = 0;
      symb = bin2dec(bits.mid(i*k, k));
      while (bitmap(pos) != symb) { pos++; }
      output(i) = symbols(pos);
    }
    return output;
  }

  ivec Modulator_1d::demodulate(const vec &signal) const
  {
    double dist, mindist;
    int closest;
    ivec output(signal.size());

    for (int i = 0; i < signal.size(); i++) {
      mindist = std::fabs(symbols(0) - signal(i));
      closest = 0;
      for (int j = 1; j < M; j++) {
	dist = std::fabs(symbols(j) - signal(i));
	if (dist < mindist) {
	  mindist = dist;
	  closest = j;
	}
      }
      output(i) = closest;
    }
    return output;
  }

  bvec Modulator_1d::demodulate_bits(const vec &signal) const
  {
    double dist, mindist;
    int closest;
    bvec output(k*signal.size());

    for (int i = 0; i < signal.size(); i++) {
      mindist = std::fabs(symbols(0) - signal(i));
      closest = 0;
      for (int j = 1; j < M; j++) {
	dist = std::fabs(symbols(j) - signal(i));
	if (dist < mindist) { 
	  mindist = dist;
	  closest = j;
	}
      }
      output.replace_mid(i*k, dec2bin(k, bitmap(closest)));
    }
    return output;
  }


  //------------- class: Modulator_2d ----------------
  Modulator_2d::Modulator_2d(const cvec &symbols, const ivec &bitmap)
  {
    set(symbols, bitmap);
  }

  void Modulator_2d::set(const cvec &in_symbols, const ivec &in_bitmap)
  {
    it_assert(in_symbols.size() == in_bitmap.size(), "Modulator_2d::set() Number of symbols and bitmap does not match");
    symbols = in_symbols;
    bitmap = in_bitmap; 
    M = bitmap.size();
    k = needed_bits(M - 1);
    calculate_softbit_matricies(bitmap);
  }

  cvec Modulator_2d::modulate(const ivec &symbolnumbers) const
  {
    cvec temp(symbolnumbers.length());
    for(int i = 0; i < symbolnumbers.length(); i++)
      temp(i) = symbols(symbolnumbers(i));
    return temp;
  }

  void Modulator_2d::modulate_bits(const bvec &bits, cvec &output) const
  {
    int no_symbols = bits.length() / k;
    int pos, symb;
    output.set_size(no_symbols);

    for (int i = 0; i < no_symbols; i++) {
      pos = 0;
      symb = bin2dec(bits.mid(i*k, k));
      while (bitmap(pos) != symb) { pos++; }
      output(i) = symbols(pos);
    }
  }

  cvec Modulator_2d::modulate_bits(const bvec &bits) const
  {
    cvec output;
    modulate_bits(bits, output);
    return output;
  }

  ivec Modulator_2d::demodulate(const cvec &signal) const
  {
    double dist, mindist;
    int closest;
    ivec output(signal.size());

    for (int i = 0; i < signal.size(); i++) {
      mindist = std::abs(symbols(0) - signal(i));
      closest = 0;
      for (int j = 1; j < M; j++) {
	dist = std::abs(symbols(j) - signal(i));
	if (dist < mindist) {
	  mindist = dist;
	  closest = j;
	}
      }
      output(i) = closest;
    }
    return output;
  }

  void Modulator_2d::demodulate_bits(const cvec &signal, bvec &bits) const
  {
    double dist, mindist;
    int closest;
    bits.set_size(k*signal.size());

    for (int i = 0; i < signal.size(); i++) {
      mindist = std::abs(symbols(0) - signal(i));
      closest = 0;
      for (int j = 1; j < M; j++) {
	dist = std::abs(symbols(j) - signal(i));
	if (dist < mindist) {
	  mindist = dist;
	  closest = j;
	}
      }
      bits.replace_mid(i*k, dec2bin(k, bitmap(closest)));
    }
  }

  bvec Modulator_2d::demodulate_bits(const cvec &signal) const
  {
    bvec bits;
    demodulate_bits(signal, bits);
    return bits;
  }

  void Modulator_2d::demodulate_soft_bits(const cvec &rx_symbols, double N0,
					  vec &soft_bits) const
  {
    // Definitions of local variables
    long no_symbols = rx_symbols.length();
    double P0, P1;
    vec soft_word(k);

    // Allocate storage space for the result vector:
    soft_bits.set_size(k*no_symbols, false);

    // For each symbol l do:
    for (int l = 0; l < no_symbols; l++) {
      // For each bit position i do:
      for (int i = 0; i < k; i++) {
	P0 = 0;
	P1 = 0;
	for (int j = 0; j < (M >> 1); j++) {
// 	  s0 = symbols( S0(i,j) );
// 	  s1 = symbols( S1(i,j) );
// 	  p_z_s0 = (1.0/(pi*N0)) * std::exp(-sqr(std::abs(z - s0)) / N0);
// 	  p_z_s1 = (1.0/(pi*N0)) * std::exp(-sqr(std::abs(z - s1)) / N0);
// 	  P0(i) += p_z_s0;
// 	  P1(i) += p_z_s1;
	  // sqr(complex) = sqr(abs(complex))
	  P0 += std::exp(-sqr(rx_symbols(l) - symbols(S0(i, j))) / N0);
	  P1 += std::exp(-sqr(rx_symbols(l) - symbols(S1(i, j))) / N0);
	}
	// The soft bits for the l-th received symbol:
	soft_word(i) = trunc_log(P0) - trunc_log(P1);
      }
      // Put the results in the result vector:
      soft_bits.replace_mid(l*k, soft_word);
    }

  }

  void Modulator_2d::demodulate_soft_bits(const cvec &rx_symbols, 
					  const cvec &chan, double N0, 
					  vec &soft_bits) const
  {
    // Definitions of local variables
    long no_symbols = rx_symbols.length();
    double a2, P0, P1;
    vec soft_word(k);

    //Allocate storage space for the result vector:
    soft_bits.set_size(k*no_symbols,false);

    // For each symbol l do:
    for (int l = 0; l < no_symbols; l++) {
      a2 = sqr(chan(l));
      // For each bit position i do:
      for (int i = 0; i < k; i++) {
	P0 = 0;
	P1 = 0;
	for (int j = 0; j < (M >> 1); j++) {
// 	  s0 = symbols( S0(i,j) );
// 	  s1 = symbols( S1(i,j) );
// 	  p_z_s0 = (1.0/(pi*N0*a2)) * std::exp(-sqr(std::abs(z-a2*s0)) / (N0*a2));
// 	  p_z_s1 = (1.0/(pi*N0*a2)) * std::exp(-sqr(std::abs(z-a2*s1)) / (N0*a2));
// 	  P0(i) += p_z_s0;
// 	  P1(i) += p_z_s1;
	  // sqr(complex) = sqr(abs(complex))
	  P0 += std::exp(-sqr(rx_symbols(l) - a2*symbols(S0(i, j))) / (N0*a2));
	  P1 += std::exp(-sqr(rx_symbols(l) - a2*symbols(S1(i, j))) / (N0*a2));
	}
	// The soft bits for the l-th received symbol:
	soft_word(i) = trunc_log(P0) - trunc_log(P1);
      }
      // Put the results in the result vector:
      soft_bits.replace_mid(l*k, soft_word);
    }
  }

  void Modulator_2d::demodulate_soft_bits_approx(const cvec &rx_symbols, 
						 double N0, 
						 vec &soft_bits) const
  {
    // Definitions of local variables
    long no_symbols = rx_symbols.length();
    double d_0, d_1;
    vec d(M);

    // Allocate storage space for the result vector:
    soft_bits.set_size(k*no_symbols, false);

    // For each symbol l do:
    for (int l = 0; l < no_symbols; l++) {
      // Calculate all distances:
      for (int i = 0; i < M; i++) { 
	d(i) = std::abs(rx_symbols(l) - symbols(i));
      }
      //For each of the k bits do:
      for (int i = 0; i < k; i++) {
	// Find the closest 0-point and the closest 1-point:
	d_0 = d(S0(i, 0));
	d_1 = d(S1(i, 0));
	for (int j = 1; j < (M >> 1); j++) {
	  if (d(S0(i, j)) < d_0) { d_0 = d(S0(i, j)); }
	  if (d(S1(i, j)) < d_1) { d_1 = d(S1(i, j)); }
	}
	//calculate the approximative metric:
	soft_bits(l*k+i) = (sqr(d_1) - sqr(d_0)) / N0;
      }
    }
  }

  void Modulator_2d::demodulate_soft_bits_approx(const cvec &rx_symbols,
						 const cvec &chan, double N0,
						 vec &soft_bits) const
  {
    // Definitions of local variables
    long no_symbols = rx_symbols.length();
    double d_0, d_1, Kf, c2;
    vec d(M);

    // Allocate storage space for the result vector:
    soft_bits.set_size(k*no_symbols, false);

    // For each symbol l do:
    for (int l = 0; l < no_symbols; l++) {
      c2 = sqr(chan(l)); 
      Kf = 1.0 / ( c2 * N0 );   
      // Calculate all distances:
      for (int i = 0; i < M; i++) {
	d(i) = std::abs(rx_symbols(l) - c2 * symbols(i));
      }
      // For each of the k bits do:
      for (int i = 0; i < k; i++) {
	//Find the closest 0-point and the closest 1-point:
	d_0 = d(S0(i, 0));
	d_1 = d(S1(i, 0));
	for (int j = 1; j < (M >> 1); j++) {
	  if (d(S0(i, j)) < d_0) { d_0 = d(S0(i, j)); }
	  if (d(S1(i, j)) < d_1) { d_1 = d(S1(i, j)); }
	}
	//calculate the approximative metric:
	soft_bits(l*k+i) = Kf * (sqr(d_1) - sqr(d_0));
      }
    }
  }

  void Modulator_2d::calculate_softbit_matricies(ivec inbitmap)
  {
    // Definitions of local variables
    int count0, count1;
    bvec bits(k);

    // Allocate storage space for the result matricies:
    S0.set_size(k, M/2, false);
    S1.set_size(k, M/2, false);

    for (int kk = 0; kk < k; kk++) {
      count0 = 0; 
      count1 = 0;
      for (int m = 0; m < M; m++) {
	bits = dec2bin(k, inbitmap(m));
	if (bits(kk) == bin(0)) {
	  S0(kk, count0) = m;
	  count0++; 
	} else {
	  S1(kk, count1) = m;
	  count1++;
	}
      }
    }
  }


  //------------- class: BPSK ----------------

  void BPSK::modulate_bits(const bvec &bits, vec &out) const
  {
    out.set_size(bits.size(), false);
    for (int i = 0; i < bits.size(); i++) { 
      out(i) = (bits(i) == 0 ? 1.0 : -1.0);
    }
  }

  // output is complex but symbols in real part only
  void BPSK::modulate_bits(const bvec &bits, cvec &out) const
  {
    out.set_size(bits.size(), false);
    for (int i = 0; i < bits.size(); i++) { 
      out(i) = (bits(i) == 0 ? 1.0 : -1.0);
    }
  }

  cvec BPSK::modulate_bits(const bvec &bits) const
  {
    cvec temp(bits.size());
    modulate_bits(bits, temp);
    return temp;
  }

  void BPSK::demodulate_bits(const vec &signal, bvec &out) const
  {
    out.set_size(signal.size(), false);
    for (int i = 0; i < signal.length(); i++) { 
      out(i) = (signal(i) > 0) ? bin(0) : bin(1); 
    }
  }

  bvec BPSK::demodulate_bits(const vec &signal) const
  {
    bvec temp(signal.size());
    demodulate_bits(signal, temp);
    return temp;
  }

  // Symbols are in real part. Channel estimation already applied to the
  // signal
  void BPSK::demodulate_bits(const cvec &signal, bvec &out) const
  {
    out.set_size(signal.size(), false);
    for (int i = 0; i < signal.length(); i++) { 
      out(i) = (std::real(signal(i)) > 0) ? bin(0) : bin(1); 
    }
  }

  bvec BPSK::demodulate_bits(const cvec &signal) const
  {
    bvec temp(signal.size());
    demodulate_bits(signal, temp);
    return temp;
  }

  // Outputs log-likelihood ratio of log(Pr(b=0|rx_symbols)/Pr(b=1|rx_symbols))
  void BPSK::demodulate_soft_bits(const vec &rx_symbols, double N0,
				  vec &soft_bits) const
  {
    double factor = 4 / N0;
    soft_bits.set_size(rx_symbols.size(), false);

    for (int i = 0; i < rx_symbols.size(); i++) {
      soft_bits(i) = factor * rx_symbols(i);
    }
  }

  // Outputs log-likelihood ratio of log(Pr(b=0|rx_symbols)/Pr(b=1|rx_symbols))
  void BPSK::demodulate_soft_bits(const cvec &rx_symbols, double N0,
				  vec &soft_bits) const
  {
    double factor = 4 / N0;
    soft_bits.set_size(rx_symbols.size(), false);

    for (int i = 0; i < rx_symbols.size(); i++) {
      soft_bits(i) = factor * std::real(rx_symbols(i));
    }
  }

  // Outputs log-likelihood ratio for fading channels
  void BPSK::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
				  double N0, vec &soft_bits) const
  {
    double factor = 4 / N0;
    soft_bits.set_size(rx_symbols.size(), false);

    for (int i = 0; i < rx_symbols.size(); i++) {
      soft_bits(i) = factor * std::real(rx_symbols(i) * std::conj(channel(i)));
    }
  }

  void BPSK::demodulate_soft_bits_approx(const cvec &rx_symbols, double N0,
					 vec &soft_bits) const
  {
    demodulate_soft_bits(rx_symbols, N0, soft_bits);
  }

  void BPSK::demodulate_soft_bits_approx(const cvec &rx_symbols, 
					 const cvec &channel, double N0,
					 vec &soft_bits) const
  {
    demodulate_soft_bits(rx_symbols, channel, N0, soft_bits);
  }



  //------------- class: PAM ----------------

  void PAM::modulate_bits(const bvec &bits, vec &out) const
  {
    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits.length() % k) {
      it_warning("PAM::modulate_bits(): The number of input bits is not a multiple of log2(M), where M is a constellation size. Remainder bits are not modulated.");
    }
    int no_symbols = bits.length() / k;

    out.set_size(no_symbols, false);
    for (int i = 0; i < no_symbols; i++) {
      out(i) = symbols(bits2symbols(bin2dec(bits.mid(i*k, k))));
    }
  }

  void PAM::modulate_bits(const bvec &bits, cvec &out) const
  {
    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits.length() % k) {
      it_warning("PAM::modulate_bits(): The number of input bits is not a multiple of log2(M), where M is a constellation size. Remainder bits are not modulated.");
    }
    int no_symbols = bits.length() / k;

    out.set_size(no_symbols, false);
    for (int i = 0; i < no_symbols; i++) {
      out(i) = symbols(bits2symbols(bin2dec(bits.mid(i*k, k))));
    }
  }

  cvec PAM::modulate_bits(const bvec &bits) const
  {
    cvec temp(bits.size());
    modulate_bits(bits, temp);
    return temp;
  }

  void PAM::demodulate_bits(const vec &signal, bvec &out) const
  {
    int est_symbol;
    out.set_size(k*signal.size(), false);

    for (int i = 0; i < signal.size(); i++) {
      est_symbol = round_i((M-1) - (signal(i) * scaling_factor + (M-1)) / 2);
      if (est_symbol < 0) 
	est_symbol = 0;
      else if (est_symbol > (M-1)) 
	est_symbol = M-1;
      out.replace_mid(i*k, bitmap.get_row(est_symbol));
    }
  }

  void PAM::demodulate_bits(const cvec &signal, bvec &out) const
  {
    int est_symbol;
    out.set_size(k*signal.size(), false);

    for (int i = 0; i < signal.size(); i++) {
      est_symbol = round_i((M-1) - (std::real(signal(i)) * scaling_factor 
				    + (M-1)) / 2);
      if (est_symbol < 0)
	est_symbol = 0;
      else if (est_symbol > (M-1))
	est_symbol = M-1;
      out.replace_mid(i*k, bitmap.get_row(est_symbol));
    }
  }

  bvec PAM::demodulate_bits(const cvec &signal) const
  {
    bvec temp(signal.size());
    demodulate_bits(signal, temp);
    return temp;
  }


  void PAM::demodulate_soft_bits(const cvec &rx_symbols, double N0,
				 vec &soft_bits) const
  {
    double P0, P1;
    vec expd(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      // calculate exponent of metric for each constellation point
      for (int j = 0; j < M; j++) {
	expd(j) = std::exp(-sqr(std::real(rx_symbols(l) - symbols(j))) / N0);
      }

      for (int i = 0; i < k; i++) {
	P0 = 0;
	P1 = 0;
	for (int j = 0; j < (M >> 1); j++) {
	  P0 += expd(S0(i, j));
	  P1 += expd(S1(i, j)); 
	}
	soft_bits(l*k+i) = trunc_log(P0) - trunc_log(P1);
      }
    }
  }

  void PAM::demodulate_soft_bits_approx(const cvec &rx_symbols, double N0,
					vec &soft_bits) const
  {
    double d0min, d1min;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      // calculate exponent of metric for each constellation point
      for (int j = 0; j < M; j++) {
	d(j) = sqr(std::real(rx_symbols(l) - symbols(j)));
      }

      for (int i = 0; i < k; i++) {
	d0min = std::numeric_limits<double>::max();
	d1min = std::numeric_limits<double>::max();
	for (int j = 0; j < (M >> 1); j++) {
	  if (d(S0(i, j)) < d0min) { d0min = d(S0(i, j)); }
	  if (d(S1(i, j)) < d1min) { d1min = d(S1(i, j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min) / N0;
      }
    }
  }

  void PAM::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
				 double N0, vec &soft_bits) const
  {
    double P0, P1;
    vec expd(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      // calculate exponent of metric for each constellation point
      for (int j = 0; j < M; j++) {
	expd(j) = std::exp(-sqr(std::real(rx_symbols(l) 
					  - channel(l) * symbols(j))) / N0);
      }

      for (int i = 0; i < k; i++) {
	P0 = 0;
	P1 = 0;
	for (int j = 0; j < (M >> 1); j++) {
	  P0 += expd(S0(i, j));
	  P1 += expd(S1(i, j)); 
	}
	soft_bits(l*k+i) = trunc_log(P0) - trunc_log(P1);
      }
    }
  }

  void PAM::demodulate_soft_bits_approx(const cvec &rx_symbols, 
					const cvec &channel, double N0,
					vec &soft_bits) const
  {
    double d0min, d1min;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      // calculate exponent of metric for each constellation point
      for (int j = 0; j < M; j++) {
	d(j) = sqr(std::real(rx_symbols(l) - channel(l) * symbols(j)));
      }

      for (int i = 0; i < k; i++) {
	d0min = std::numeric_limits<double>::max();
	d1min = std::numeric_limits<double>::max();
	for (int j = 0; j < (M >> 1); j++) {
	  if (d(S0(i, j)) < d0min) { d0min = d(S0(i, j)); }
	  if (d(S1(i, j)) < d1min) { d1min = d(S1(i, j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min) / N0;
      }
    }
  }

  void PAM::set_M(int Mary)
  {
    M = Mary;
    k = needed_bits(M - 1);

    it_assert(pow2i(k) == Mary, "PAM::set_M(): M is not a power of 2");

    int count0, count1;
    bvec bits;

    symbols.set_size(M, false);
    bits2symbols.set_size(M, false);
    bitmap = graycode(k);
    average_energy = (sqr(M) - 1) / 3.0;
    scaling_factor = std::sqrt(average_energy);

    for (int i = 0; i < M; i++) {
      symbols(i) = ((M-1) - i*2) / scaling_factor;
      bits2symbols(bin2dec(bitmap.get_row(i))) = i;
    }

    // Calculate the soft bit mapping matrices S0 and S1
    S0.set_size(k, M/2, false);
    S1.set_size(k, M/2, false);
  
    for (int kk = 0; kk < k; kk++) {
      count0 = 0; 
      count1 = 0;
      for (int m = 0; m < M; m++) {
	bits = bitmap.get_row(m);
	if (bits(kk) == bin(0)) {
	  S0(kk, count0) = m;
	  count0++; 
	} else {
	  S1(kk, count1) = m;
	  count1++;
	}
      }
    }
  }


  //------------- class: QPSK ----------------

  void QPSK::modulate_bits(const bvec &bits, cvec &out) const
  {
    // Check if some bits have to be cut and print warning message in such
    // case.
    if (bits.length() & 0x1) {
      it_warning("QPSK::modulate_bits(): The number of input bits is not a multiple of 2. Remainder bits are not modulated.");
    }
    double real_part, imag_part;
    int no_symbols = bits.length() >> 1;

    out.set_size(no_symbols, false);
    for (int i = 0; i < no_symbols; i++) {
      real_part = (bits(2*i) == 0) ? M_SQRT1_2 : -M_SQRT1_2;
      imag_part = (bits(2*i+1) == 0) ? M_SQRT1_2 : -M_SQRT1_2;
      out(i) = std::complex<double>(real_part, imag_part);
    }
  }

  cvec QPSK::modulate_bits(const bvec &bits) const
  {
    cvec out;
    modulate_bits(bits, out);
    return out;
  }

  void QPSK::demodulate_bits(const cvec &signal, bvec &out) const
  {
    int no_symbols = signal.size();
    out.set_size(2*no_symbols, false);
    for (int i = 0; i < no_symbols; i++) {
      out(2*i) = ((std::real(signal(i)) > 0) ? 0 : 1);
      out(2*i+1) = ((std::imag(signal(i)) > 0) ? 0 : 1);
    }
  }

  bvec QPSK::demodulate_bits(const cvec &signal) const
  {
    bvec out;
    demodulate_bits(signal, out);
    return out;
  }


  // Outputs log-likelihood ratio of log (Pr(b=0|rx_symbols)/Pr(b=1|rx_symbols))
  void QPSK::demodulate_soft_bits(const cvec &rx_symbols, double N0,
				  vec &soft_bits) const
  {
    soft_bits.set_size(2*rx_symbols.size(), false);
    double factor = 2 * std::sqrt(2.0) / N0;
    for (int i = 0; i < rx_symbols.size(); i++) {
      soft_bits(2*i) = std::real(rx_symbols(i)) * factor;
      soft_bits(2*i+1) = std::imag(rx_symbols(i)) * factor;
    }
  }

  // Outputs log-likelihood ratio for fading channels
  void QPSK::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
				  double N0, vec &soft_bits) const
  {
    soft_bits.set_size(2*rx_symbols.size(), false);
    std::complex<double> temp;
    double factor = 2 * std::sqrt(2.0) / N0;
    
    for (int i = 0; i < rx_symbols.size(); i++) {
      temp = rx_symbols(i) * std::conj(channel(i));
      soft_bits(2*i) = std::real(temp) * factor;
      soft_bits(2*i+1) = std::imag(temp) * factor;
    }
  }


  void QPSK::demodulate_soft_bits_approx(const cvec &rx_symbols, double N0,
					 vec &soft_bits) const
  {
    demodulate_soft_bits(rx_symbols, N0, soft_bits);
  }

  void QPSK::demodulate_soft_bits_approx(const cvec &rx_symbols,
					 const cvec &channel, double N0,
					 vec &soft_bits) const
  {
    demodulate_soft_bits(rx_symbols, channel, N0, soft_bits);
  }


  //------------- class: PSK ----------------

  void PSK::modulate_bits(const bvec &bits, cvec &out) const
  {
    // Check if some bits have to be cut and print warning message in such
    // case.
    if (bits.length() % k) {
      it_warning("PSK::modulate_bits(): The number of input bits is not a multiple of log2(M), where M is a constellation size. Remainder bits are not modulated.");
    }
    int no_symbols = bits.length() / k;

    out.set_size(no_symbols, false);

    for (int i = 0; i < no_symbols; i++) {
      out(i) = symbols(bits2symbols(bin2dec(bits.mid(i*k, k))));
    }
  }

  cvec PSK::modulate_bits(const bvec &bits) const
  {
    cvec temp(bits.size());
    modulate_bits(bits, temp);
    return temp;
  }

  void PSK::demodulate_bits(const cvec &signal, bvec &out) const
  {
    int est_symbol;
    double ang, temp;

    out.set_size(k*signal.size(), false);

    for (int i = 0; i < signal.size(); i++) {
      ang = std::arg(signal(i));
      temp = (ang < 0) ? (2*pi+ang) : ang;
      est_symbol = round_i(temp * (M/2) / pi) % M;
      out.replace_mid(i*k, bitmap.get_row(est_symbol));
    }
  }

  bvec PSK::demodulate_bits(const cvec &signal) const
  {
    bvec temp;
    demodulate_bits(signal, temp);
    return temp;
  }

  void PSK::demodulate_soft_bits(const cvec &rx_symbols, double N0,
				 vec &soft_bits) const
  {
    double P0, P1;
    vec expd(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
	expd(j) = std::exp(-sqr(rx_symbols(l) - symbols(j)) / N0);
      }
      for (int i = 0; i < k; i++) {
	P0 = 0;
	P1 = 0;
	for (int j = 0; j < (M >> 1); j++) {
	  P0 += expd(S0(i, j));
	  P1 += expd(S1(i, j));  
	}
	soft_bits(l*k+i) = trunc_log(P0) - trunc_log(P1);
      }
    }
  }

  void PSK::demodulate_soft_bits_approx(const cvec &rx_symbols, double N0,
					vec &soft_bits) const
  {
    double d0min, d1min;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++)
	d(j) = sqr(rx_symbols(l) - symbols(j));
      for (int i = 0; i < k; i++) {
	d0min = std::numeric_limits<double>::max();
	d1min = std::numeric_limits<double>::max();
	for (int j = 0; j < (M >> 1); j++) {
	  if (d(S0(i, j)) < d0min) { d0min = d(S0(i, j)); }
	  if (d(S1(i, j)) < d1min) { d1min = d(S1(i, j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min) / N0;
      }
    }
  }

  void PSK::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
				 double N0, vec &soft_bits) const
  {
    double P0, P1;
    vec expd(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
	expd(j) = std::exp(-sqr(rx_symbols(l) - channel(l) * symbols(j)) / N0);
      }
      for (int i = 0; i < k; i++) {
	P0 = 0;
	P1 = 0;
	for (int j = 0; j < (M >> 1); j++) {
	  P0 += expd(S0(i, j));
	  P1 += expd(S1(i, j));  
	}
	soft_bits(l*k+i) = trunc_log(P0) - trunc_log(P1);
      }
    }
  }

  void PSK::demodulate_soft_bits_approx(const cvec &rx_symbols,
					const cvec &channel, double N0,
					vec &soft_bits) const
  {
    double d0min, d1min;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++)
	d(j) = sqr(rx_symbols(l) - channel(l) * symbols(j));
      for (int i = 0; i < k; i++) {
	d0min = std::numeric_limits<double>::max();
	d1min = std::numeric_limits<double>::max();
	for (int j = 0; j < (M >> 1); j++) {
	  if (d(S0(i, j)) < d0min) { d0min = d(S0(i, j)); }
	  if (d(S1(i, j)) < d1min) { d1min = d(S1(i, j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min) / N0;
      }
    }
  }

  void PSK::set_M(int Mary)
  {
    k = needed_bits(Mary - 1);
    M = Mary;
    it_assert(pow2i(k) == Mary, "PSK::set_M(): M is not a power of 2");
    symbols.set_size(M, false);
    bits2symbols.set_size(M, false);
    bitmap = graycode(k);

    double delta = 2.0 * pi / M;
    double epsilon = delta / 10000.0;
    std::complex<double> symb;
    for (int i = 0; i < M; i++) {
      symb = std::complex<double>(std::polar(1.0, delta*i));
      if (std::fabs(std::real(symb)) < epsilon) { 
	symbols(i) = std::complex<double>(0.0, std::imag(symb)); 
      }
      else if (std::fabs(std::imag(symb)) < epsilon) { 
	symbols(i) = std::complex<double>(std::real(symb), 0.0);
      }
      else { symbols(i) = symb; }

      bits2symbols(bin2dec(bitmap.get_row(i))) = i;
    }

    // Calculate the soft bit mapping matrices S0 and S1
    S0.set_size(k, M/2,false);
    S1.set_size(k, M/2,false);
    int count0, count1;
    bvec bits;
  
    for (int kk = 0; kk < k; kk++) {
      count0 = 0; 
      count1 = 0;
      for (int m = 0; m < M; m++) {
	bits = bitmap.get_row(m);
	if (bits(kk) == bin(0)) {
	  S0(kk, count0) = m;
	  count0++; 
	} else {
	  S1(kk, count1) = m;
	  count1++;
	}
      }
    }
  }


  //------------- class: QAM ----------------

  void QAM::modulate_bits(const bvec &bits, cvec &out) const
  {
    // Check if some bits have to be cut and print warning message in such
    // case.
    if (bits.length() % k) {
      it_warning("QAM::modulate_bits(): The number of input bits is not a multiple of log2(M), where M is a constellation size. Remainder bits are not modulated.");
    }
    int no_symbols = bits.length() / k;
    out.set_size(no_symbols, false);
    for (int i = 0; i < no_symbols; i++) {
      out(i) = symbols(bits2symbols(bin2dec(bits.mid(i*k, k))));
    }
  }

  cvec QAM::modulate_bits(const bvec &bits) const
  {
    cvec temp(bits.size());
    modulate_bits(bits, temp);
    return temp;
  }

  void QAM::demodulate_bits(const cvec &signal, bvec &out) const
  {
    out.set_size(k*signal.size(), false);

    int temp_real, temp_imag;

    for (int i = 0; i < signal.size(); i++) {
      temp_real = round_i((L-1) - (std::real(signal(i) * scaling_factor)
				   + (L-1)) / 2.0);
      temp_imag = round_i((L-1) - (std::imag(signal(i) * scaling_factor)
				   + (L-1)) / 2.0);
      if (temp_real < 0) 
	temp_real = 0; 
      else if (temp_real > (L-1)) 
	temp_real = (L-1);
      if (temp_imag < 0) 
	temp_imag = 0; 
      else if (temp_imag > (L-1)) 
	temp_imag = (L-1);
      out.replace_mid(k*i, bitmap.get_row(temp_imag * L + temp_real));
    }
  }

  bvec QAM::demodulate_bits(const cvec &signal) const
  {
    bvec temp;
    demodulate_bits(signal, temp);
    return temp;
  }

  void QAM::demodulate_soft_bits(const cvec &rx_symbols, double N0,
				 vec &soft_bits) const
  {
    double P0, P1;
    vec expd(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
	expd(j) = std::exp(-sqr(rx_symbols(l) - symbols(j)) / N0);
      }
      for (int i = 0; i < k; i++) {
	P0 = 0;
	P1 = 0;
	for (int j = 0; j < (M >> 1); j++) {
	  P0 += expd(S0(i, j));
	  P1 += expd(S1(i, j));  
	}
	soft_bits(l*k+i) = trunc_log(P0) - trunc_log(P1);
      }
    }
  }


  void QAM::demodulate_soft_bits_approx(const cvec &rx_symbols, double N0,
					vec &soft_bits) const
  {
    double d0min, d1min;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++)
	d(j) = sqr(rx_symbols(l) - symbols(j));

      for (int i = 0; i < k; i++) {
	d0min = std::numeric_limits<double>::max();
	d1min = std::numeric_limits<double>::max();
	for (int j = 0; j < (M >> 1); j++) {
	  if (d(S0(i, j)) < d0min) { d0min = d(S0(i, j)); }
	  if (d(S1(i, j)) < d1min) { d1min = d(S1(i, j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min) / N0;
      }
    }
  }

  void QAM::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
				 double N0, vec &soft_bits) const
  {
    double P0, P1;
    vec expd(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++) {
	expd(j) = std::exp(-sqr(rx_symbols(l) - channel(l) * symbols(j)) / N0);
      }
      for (int i = 0; i < k; i++) {
	P0 = 0;
	P1 = 0;
	for (int j = 0; j < (M >> 1); j++) {
	  P0 += expd(S0(i, j));
	  P1 += expd(S1(i, j));  
	}
	soft_bits(l*k+i) = trunc_log(P0) - trunc_log(P1);
      }
    }
  }


  void QAM::demodulate_soft_bits_approx(const cvec &rx_symbols, 
					const cvec &channel, double N0,
					vec &soft_bits) const
  {
    double d0min, d1min;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (int l = 0; l < rx_symbols.size(); l++) {
      for (int j = 0; j < M; j++)
	d(j) = sqr(rx_symbols(l) - channel(l) * symbols(j));

      for (int i = 0; i < k; i++) {
	d0min = std::numeric_limits<double>::max();
	d1min = std::numeric_limits<double>::max();
	for (int j = 0; j < (M >> 1); j++) {
	  if (d(S0(i, j)) < d0min) { d0min = d(S0(i, j)); }
	  if (d(S1(i, j)) < d1min) { d1min = d(S1(i, j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min) / N0;
      }
    }
  }


  void QAM::set_M(int Mary)
  {
    k = needed_bits(Mary - 1);
    M = Mary;
    L = round_i(std::sqrt(static_cast<double>(M)));
    it_assert(pow2i(k) == Mary, "QAM::set_M(): M is not a power of 2");

    int count0, count1;
    bvec bits;

    symbols.set_size(M, false);
    bitmap.set_size(M, k, false);
    bits2symbols.set_size(M, false);
    bmat gray_code = graycode(needed_bits(L - 1));
    average_energy = (M-1) * 2.0 / 3.0;
    scaling_factor = std::sqrt(average_energy);

    for (int i = 0; i < L; i++) {
      for (int j = 0; j < L; j++) {
	symbols(i*L+j) = std::complex<double>(((L-1) - j*2) / scaling_factor,
					      ((L-1) - i*2) / scaling_factor);
	bitmap.set_row(i*L+j, concat(gray_code.get_row(i), 
				     gray_code.get_row(j)));
	bits2symbols(bin2dec(bitmap.get_row(i*L+j))) = i*L+j;
      }
    }

    // Calculate the soft bit mapping matrices S0 and S1
    S0.set_size(k, M/2, false);
    S1.set_size(k, M/2, false);
  
    for (int kk = 0; kk < k; kk++) {
      count0 = 0; 
      count1 = 0;
      for (int m = 0; m < M; m++) {
	bits = bitmap.get_row(m);
	if (bits(kk) == bin(0)) {
	  S0(kk, count0) = m;
	  count0++; 
	} else {
	  S1(kk, count1) = m;
	  count1++;
	}
      }
    }
  }

} // namespace itpp
