/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2005 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
 * Tobias Ringström, and Jonas Samuelsson.                                   *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*! 
  \file 
  \brief Implementation of modulator classes
  \author 

  1.19

  2003/06/17 11:48:36
*/

#include <float.h>

#include "itpp/base/binary.h"
#include "itpp/base/matfunc.h"
#include "itpp/base/specmat.h"
#include "itpp/comm/modulator.h"
#include "itpp/comm/commfunc.h"

#include <complex>


namespace itpp {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define M_1_DIV_SQRT2 0.70710678118654752440
#endif //DOXYGEN_SHOULD_SKIP_THIS

  //------------- class: Modulator_1d ----------------

  Modulator_1d::Modulator_1d(const vec &insymbols, const ivec &inbitmap)
  {
    set(insymbols, inbitmap);
  }

  void Modulator_1d::set(const vec &insymbols, const ivec &inbitmap)
  {
    it_assert(insymbols.size() == inbitmap.size(), "Modulator_1d: number of symbols and bitmap does not match");
    M = inbitmap.size();
    k = round_i(log2(double(M)));
    symbols = insymbols;
    bitmap = inbitmap; 
  }

  vec Modulator_1d::modulate(const ivec &symbolnumbers) const
  {
    vec temp(symbolnumbers.size());
    for (int i=0; i<symbolnumbers.size(); i++)
      temp(i) = symbols(symbolnumbers(i));
    return temp;
  }

  vec Modulator_1d::modulate_bits(const bvec &bits) const
  {
    int no_symbols = floor_i(double(bits.length())/double(k)), i, pos, symb;
    vec output = zeros(no_symbols);

    for (i=0; i<no_symbols; i++) {
      pos = 0;
      symb = bin2dec(bits.mid(i*k,k));
      while (bitmap(pos)!=symb) { pos++; }
      output(i) = symbols(pos);
    }
    return output;
  }

  ivec Modulator_1d::demodulate(const vec &signal) const
  {
    int i, j;
    double dist, mindist;
    int closest;
    ivec output(signal.size());

    for (i=0; i<signal.size(); i++) {
      mindist = std::abs(double(symbols(0)) - signal(i));
      closest = 0;
      for (j=1; j<M; j++) {
	dist = std::abs(double(symbols(j)) - signal(i));
	if (dist<mindist) { mindist = dist; closest = j; }
      }
      output(i) = closest;
    }
    return output;
  }

  bvec Modulator_1d::demodulate_bits(const vec &signal) const
  {
    int i, j;
    double dist, mindist;
    int closest;
    bvec output(k*signal.size());

    for (i=0; i<signal.size(); i++) {
      mindist = std::abs(double(symbols(0)) - signal(i));
      closest = 0;
      for (j=1; j<M; j++) {
	dist = std::abs(double(symbols(j)) - signal(i));
	if (dist<mindist){ mindist = dist; closest = j; }
      }
      output.replace_mid(i*k,dec2bin(k,bitmap(closest)));
    }
    return output;
  }

  vec Modulator_1d::get_symbols(void) const { return symbols; }

  ivec Modulator_1d::get_bitmap(void) const { return bitmap; }

  //------------- class: BPSK ----------------

  void BPSK::modulate_bits(const bvec &bits, vec &out) const
  {
    out.set_size(bits.size(), false);
    for (int i=0;i<bits.size();i++) { out(i) = bits(i) == 0 ? 1.0 : -1.0;}
  }

  // output is complex but symbols in real part only
  void BPSK::modulate_bits(const bvec &bits, cvec &out) const
  {
    out.set_size(bits.size(), false);
    for (int i=0;i<bits.size();i++) { out(i) = bits(i) == 0 ? 1.0 : -1.0;}
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
    for (int i=0; i<signal.length(); i++) { out(i) = (signal(i)>0) ? bin(0) : bin(1); }
  }

  bvec BPSK::demodulate_bits(const vec &signal) const
  {
    bvec temp(signal.size());
    demodulate_bits(signal, temp);
    return temp;
  }

  // Symbols are in real part. That is channel estimation is already applied to the signal
  void BPSK::demodulate_bits(const cvec &signal, bvec &out) const
  {
    out.set_size(signal.size(), false);
    for (int i=0; i<signal.length(); i++) { out(i) = (std::real(signal(i))>0) ? bin(0) : bin(1); }
  }

  bvec BPSK::demodulate_bits(const cvec &signal) const
  {
    bvec temp(signal.size());
    demodulate_bits(signal, temp);
    return temp;
  }

  // Outputs log-likelihood ratio of log (Pr(bit = 0|rx_symbols)/Pr(bit = 1|rx_symbols))
  void BPSK::demodulate_soft_bits(const vec &rx_symbols, const double N0, vec &soft_bits) const
  {
    double factor = 4/N0;
    soft_bits.set_size(rx_symbols.size(), false);

    for (int i=0; i<rx_symbols.size(); i++) {
      soft_bits(i) = factor*rx_symbols(i);
    }
  }

  // Outputs log-likelihood ratio of log (Pr(bit = 0|rx_symbols)/Pr(bit = 1|rx_symbols))
  void BPSK::demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits) const
  {
    double factor = 4/N0;
    soft_bits.set_size(rx_symbols.size(), false);

    for (int i=0; i<rx_symbols.size(); i++) {
      soft_bits(i) = factor*std::real(rx_symbols(i));
    }
  }

  // Outputs log-likelihood ratio for fading channels
  void BPSK::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const
  {
    double factor = 4/N0;
    soft_bits.set_size(rx_symbols.size(), false);

    for (int i=0; i<rx_symbols.size(); i++) {
      soft_bits(i) = factor*std::real(rx_symbols(i)*std::conj(channel(i)));
    }
  }

  void BPSK::demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits) const
  {
    demodulate_soft_bits(rx_symbols, N0, soft_bits);
  }

  void BPSK::demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const
  {
    demodulate_soft_bits(rx_symbols, channel, N0, soft_bits);
  }


  //------------- class: PAM ----------------

  void PAM::modulate_bits(const bvec &bits, vec &out) const
  {
    int no_symbols, i, symb;
    int bits_len = bits.length();

    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits_len % k) {
      it_warning("PAM::modulate_bits(): The number of input bits is not a multiple of log2(M),\nwhere M is a constellation size. Remainder bits are not modulated.");
    }
    no_symbols = floor_i(double(bits_len) / double(k));

    out.set_size(no_symbols, false);

    for (i=0; i<no_symbols; i++) {
      symb = bin2dec(bits.mid(i*k,k));
      out(i) = symbols(bits2symbols(symb));
    }
  }

  void PAM::modulate_bits(const bvec &bits, cvec &out) const
  {
    int no_symbols, i, symb;
    int bits_len = bits.length();

    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits_len % k) {
      it_warning("PAM::modulate_bits(): The number of input bits is not a multiple of log2(M),\nwhere M is a constellation size. Remainder bits are not modulated.");
    }
    no_symbols = floor_i(double(bits_len) / double(k));

    out.set_size(no_symbols, false);

    for (i=0; i<no_symbols; i++) {
      symb = bin2dec(bits.mid(i*k,k));
      out(i) = symbols(bits2symbols(symb));
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
    int i, est_symbol;
    out.set_size(k*signal.size(), false);

    for (i=0; i<signal.size(); i++) {
      est_symbol = round_i((M-1)-(signal(i)*scaling_factor+(M-1))/2.0);
      if (est_symbol<0) est_symbol = 0;
      else if (est_symbol>(M-1)) est_symbol = (M-1);
      out.replace_mid(i*k,bitmap.get_row(est_symbol));
    }
  }

  void PAM::demodulate_bits(const cvec &signal, bvec &out) const
  {
    int i, est_symbol;
    out.set_size(k*signal.size(), false);

    for (i=0; i<signal.size(); i++) {
      est_symbol = round_i((M-1)-(std::real(signal(i))*scaling_factor+(M-1))/2.0);
      if (est_symbol<0) est_symbol = 0;
      else if (est_symbol>(M-1)) est_symbol = (M-1);
      out.replace_mid(i*k,bitmap.get_row(est_symbol));
    }
  }

  bvec PAM::demodulate_bits(const cvec &signal) const
  {
    bvec temp(signal.size());
    demodulate_bits(signal, temp);
    return temp;
  }


  void PAM::demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double P0, P1;
    double d0min, d0min2, d1min, d1min2;
    double treshhold = -log(eps)*N0; // To be sure that any precision is left in the calculatation
    double inv_N0 = 1/N0;
    vec d(M), expd(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {

      for (j=0; j<M; j++) {
	d(j) = sqr(std::real(rx_symbols(l)-symbols(j)));
	expd(j) = exp((-d(j))*inv_N0);
      }

      for (i=0; i<k; i++) {
	P0 = P1 = 0;
	d0min = d0min2 = d1min = d1min2 = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min2 = d0min; d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min2 = d1min; d1min = d(S1(i,j)); }
	  P0 += expd(S0(i,j));
	  P1 += expd(S1(i,j)); 
	}
	if ( (d0min2-d0min) > treshhold && (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
	else if ( (d0min2-d0min) > treshhold )
	  soft_bits(l*k+i) = -d0min*inv_N0-log(P1);
	else if ( (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = log(P0) + d1min*inv_N0;
	else
	  soft_bits(l*k+i) = log(P0/P1);
      }
    }
  }

  void PAM::demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double d0min, d1min;
    double inv_N0 = 1/N0;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {

      for (j=0; j<M; j++)
	d(j) = sqr(std::real(rx_symbols(l)-symbols(j)));

      for (i=0; i<k; i++) {
	d0min = d1min = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min = d(S1(i,j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
      }
    }
  }

  void PAM::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double P0, P1;
    double d0min, d0min2, d1min, d1min2;
    double treshhold = -log(eps)*N0; // To be sure that any precision is left in the calculatation
    double inv_N0 = 1/N0;
    vec d(M), expd(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {
      for (j=0; j<M; j++) {
	d(j) = sqr(std::real(rx_symbols(l)-channel(l)*symbols(j)));
	expd(j) = exp((-d(j))*inv_N0);
      }

      for (i=0; i<k; i++) {
	P0 = P1 = 0;
	d0min = d0min2 = d1min = d1min2 = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min2 = d0min; d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min2 = d1min; d1min = d(S1(i,j)); }
	  P0 += expd(S0(i,j));
	  P1 += expd(S1(i,j)); 
	}
	if ( (d0min2-d0min) > treshhold && (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
	else if ( (d0min2-d0min) > treshhold )
	  soft_bits(l*k+i) = -d0min*inv_N0-log(P1);
	else if ( (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = log(P0) + d1min*inv_N0;
	else
	  soft_bits(l*k+i) = log(P0/P1);
      }
    }
  }

  void PAM::demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double d0min, d1min;
    double inv_N0 = 1/N0;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {

      for (j=0; j<M; j++)
	d(j) = sqr(std::real(rx_symbols(l)-channel(l)*symbols(j)));

      for (i=0; i<k; i++) {
	d0min = d1min = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min = d(S1(i,j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
      }
    }
  }

  void PAM::set_M(int Mary)
  {
    M = Mary;
    k = round_i(log2(double(M)));

    it_error_if(abs(pow2i(k)-Mary)>0.0,"M-ary PAM: M is not a power of 2");

    int i, kk, m, count0, count1;
    bvec bits;

    symbols.set_size(M, false);
    bits2symbols.set_size(M, false);
    bitmap = graycode(k);
    average_energy = double(M*M-1)/3.0;
    scaling_factor = sqrt(average_energy);

    for (i=0; i<M; i++) {
      symbols(i) = double((M-1)-i*2) / scaling_factor;
      bits2symbols(bin2dec(bitmap.get_row(i))) = i;
    }

    //Calculate the soft bit mapping matrices S0 and S1
    S0.set_size(k,M/2,false);
    S1.set_size(k,M/2,false);
  
    for (kk=0; kk<k; kk++) {
      count0 = 0; 
      count1 = 0;
      for (m=0; m<M; m++) {
	bits = bitmap.get_row(m);
	if (bits(kk)==bin(0)) {
	  S0(kk,count0) = m;
	  count0++; 
	} else {
	  S1(kk,count1) = m;
	  count1++;
	}
      }
    }
  }

  //------------- class: Modulator_2d ----------------
  Modulator_2d::Modulator_2d(const cvec &insymbols, const ivec &inbitmap)
  {
    set(insymbols, inbitmap);
  }

  void Modulator_2d::set(const cvec &insymbols, const ivec &inbitmap)
  {
    it_assert(insymbols.size() == inbitmap.size(), "Modulator_2d: number of symbols and bitmap does not match");
    symbols = insymbols;
    bitmap = inbitmap; 
    M = bitmap.size();
    k = round_i(log2(double(M)));
    soft_bit_mapping_matrices_calculated = false;
  }

  cvec Modulator_2d::modulate(const ivec &symbolnumbers) const
  {
    cvec temp(symbolnumbers.length());
    for(int i=0;i<symbolnumbers.length();i++)
      temp(i)=symbols(symbolnumbers(i));
    return temp;
  }

  cvec Modulator_2d::modulate_bits(const bvec &bits) const
  {
    int no_symbols = floor_i(double(bits.length())/double(k)), i, pos, symb;
    cvec output = zeros_c(no_symbols);

    for (i=0; i<no_symbols; i++) {
      pos = 0;
      symb = bin2dec(bits.mid(i*k,k));
      while (bitmap(pos)!=symb) { pos++; }
      output(i) = symbols(pos);
    }
    return output;
  }

  ivec Modulator_2d::demodulate(const cvec &signal) const
  {
    int i, j;
    double dist, mindist;
    int closest;
    ivec output(signal.size());

    for (i=0; i<signal.size(); i++) {
      mindist = std::abs( symbols(0) - signal(i) );
      closest = 0;
      for (j=1; j<M; j++) {
	dist = std::abs( symbols(j) - signal(i) );
	if (dist<mindist){ mindist = dist; closest = j; }
      }
      output(i) = closest;
    }
    return output;
  }

  bvec Modulator_2d::demodulate_bits(const cvec &signal) const
  {
    int i, j;
    double dist, mindist;
    int closest;
    bvec output(k*signal.size());

    for (i=0; i<signal.size(); i++) {
      mindist = std::abs( symbols(0) - signal(i) );
      closest = 0;
      for (j=1; j<M; j++) {
	dist = std::abs( symbols(j) - signal(i) );
	if (dist<mindist){ mindist = dist; closest = j; }
      }
      output.replace_mid(i*k,dec2bin(k,bitmap(closest))); 
    }
    return output;
  }

  cvec Modulator_2d::get_symbols() const { return symbols; }

  ivec Modulator_2d::get_bitmap(void) const { return bitmap; }

  void Modulator_2d::demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits)
  {
    //Definitions of local variables
    long no_symbols = rx_symbols.length();
    long l, i, j;
    double p_z_s0, p_z_s1;
    std::complex<double> z, s0, s1;
    vec soft_word(k), P0(k), P1(k);

    //Check if the soft bit mapping matrices S0 and S1 needs to be calculated
    if (soft_bit_mapping_matrices_calculated==false) {
      calculate_softbit_matricies(bitmap);
    }

    //Allocate storage space for the result vector:
    soft_bits.set_size(k*no_symbols,false);

    //For each symbol l do:
    for (l=0; l<no_symbols; l++) {

      P0.clear();
      P1.clear();
      z = rx_symbols(l);

      //For each bit position i do:
      for (i=0; i<k; i++) {

	for (j=0; j<(M/2); j++) {
	  s0 = symbols( S0(i,j) );
	  s1 = symbols( S1(i,j) );
	  p_z_s0 = (1.0/(pi*N0)) * exp(-(pow(std::abs(z-s0),2.0))/N0);
	  p_z_s1 = (1.0/(pi*N0)) * exp(-(pow(std::abs(z-s1),2.0))/N0);
	  P0(i) += p_z_s0;
	  P1(i) += p_z_s1;
	}
	//The soft bits for the l-th received symbol:
	soft_word(i) = log( P0(i) / P1(i) );
      }
      //Put the results in the result vector:
      soft_bits.replace_mid(l*k,soft_word);
    }

  }

  void Modulator_2d::demodulate_soft_bits(const cvec &rx_symbols, const cvec &chan, const double N0, vec &soft_bits)
  {
    //Definitions of local variables
    long no_symbols = rx_symbols.length();
    long l, i, j;
    double p_z_s0, p_z_s1, a2;
    std::complex<double> z, s0, s1;
    vec soft_word(k), P0(k), P1(k);

    //Check if the soft bit mapping matrices S0 and S1 needs to be calculated
    if (soft_bit_mapping_matrices_calculated==false) {
      calculate_softbit_matricies(bitmap);
    }

    //Allocate storage space for the result vector:
    soft_bits.set_size(k*no_symbols,false);

    //For each symbol do:
    for (l=0; l<no_symbols; l++) {

      P0.clear();
      P1.clear();
      z = rx_symbols(l);
      a2 = pow(std::abs(chan(l)),2.0);

      //For each bit position i do:
      for (i=0; i<k; i++) {

	for (j=0; j<(M/2); j++) {
	  s0 = symbols( S0(i,j) );
	  s1 = symbols( S1(i,j) );
	  p_z_s0 = (1.0/(pi*N0*a2)) * exp(-(pow(std::abs(z-a2*s0),2.0))/(N0*a2));
	  p_z_s1 = (1.0/(pi*N0*a2)) * exp(-(pow(std::abs(z-a2*s1),2.0))/(N0*a2));
	  P0(i) += p_z_s0;
	  P1(i) += p_z_s1;
	}
	//The soft bits for the l-th received symbol:
	soft_word(i) = log( P0(i) / P1(i) );
      }
      //Put the results in the result vector:
      soft_bits.replace_mid(l*k,soft_word);
    }

  }

  void Modulator_2d::demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits)
  {
    //Definitions of local variables
    long no_symbols = rx_symbols.length();
    long l, i, j;
    double d_0, d_1, Kf;
    vec d(M);
    Kf = (1.0/N0);

    //Check if the soft bit mapping matrices S0 and S1 needs to be calculated
    if (soft_bit_mapping_matrices_calculated==false) {
      calculate_softbit_matricies(bitmap);
    }

    //Allocate storage space for the result vector:
    soft_bits.set_size(k*no_symbols,false);

    //for each symbol l do:
    for (l=0; l<no_symbols; l++) {

      //Calculate all distances:
      for (i=0; i<M; i++) { d(i) = std::abs( rx_symbols(l) - symbols(i) ); }

      //For each of the k bits do:
      for (i=0; i<k; i++) {

	//Find the closest 0-point and the closest 1-point:
	d_0 = d( S0(i,0) );
	d_1 = d( S1(i,0) );
	for (j=1; j<(M/2); j++) {
	  if ( d( S0(i,j) ) < d_0) { d_0 = d( S0(i,j) ); }
	  if ( d( S1(i,j) ) < d_1) { d_1 = d( S1(i,j) ); }
	}

	//calculate the approximative metric:
	soft_bits(l*k+i) = Kf * ( pow(d_1,2.0) - pow(d_0,2.0) );

      }

    }

  }

  void Modulator_2d::demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &chan, const double N0, vec &soft_bits)
  {

    //Definitions of local variables
    long no_symbols = rx_symbols.length();
    long l, i, j;
    double d_0, d_1, Kf, c2;
    vec d(M);

    //Check if the soft bit mapping matrices S0 and S1 needs to be calculated
    if (soft_bit_mapping_matrices_calculated==false) {
      calculate_softbit_matricies(bitmap);
    }

    //Allocate storage space for the result vector:
    soft_bits.set_size(k*no_symbols,false);

    //for each symbol l do:
    for (l=0; l<no_symbols; l++) {

      c2 = pow(std::abs(chan(l)),2.0); 
      Kf = 1.0 / ( c2 * N0 );   

      //Calculate all distances:
      for (i=0; i<M; i++) { d(i) = std::abs( rx_symbols(l) - c2*symbols(i) ); }

      //For each of the k bits do:
      for (i=0; i<k; i++) {

	//Find the closest 0-point and the closest 1-point:
	d_0 = d( S0(i,0) );
	d_1 = d( S1(i,0) );
	for (j=1; j<(M/2); j++) {
	  if ( d( S0(i,j) ) < d_0) { d_0 = d( S0(i,j) ); }
	  if ( d( S1(i,j) ) < d_1) { d_1 = d( S1(i,j) ); }
	}

	//calculate the approximative metric:
	soft_bits(l*k+i) = Kf * ( pow(d_1,2.0) - pow(d_0,2.0) );

      }

    }

  }

  void Modulator_2d::calculate_softbit_matricies(ivec inbitmap)
  {
    //Definitions of local variables
    int kk, m, count0, count1;
    bvec bits(k);

    //Allocate storage space for the result matricies:
    S0.set_size(k,M/2,false);
    S1.set_size(k,M/2,false);

    for (kk=0; kk<k; kk++) {
      count0 = 0; 
      count1 = 0;
      for (m=0; m<M; m++) {
	bits = dec2bin(k,inbitmap(m));
	if (bits(kk)==bin(0)) {
	  S0(kk,count0) = m;
	  count0++; 
	} else {
	  S1(kk,count1) = m;
	  count1++;
	}
      }
    }

  }


  //------------- class: QPSK ----------------

  void QPSK::modulate_bits(const bvec &bits, cvec &out) const
  {
    int no_symbols, i;
    int bits_len = bits.length();
    double real_part, imag_part;

    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits_len % 2) {
      it_warning("QPSK::modulate_bits(): The number of input bits is not a multiple of 2.\nRemainder bits are not modulated.");
    }
    no_symbols = floor_i(double(bits_len) / 2);

    out.set_size(no_symbols, false);
    for (i=0; i<no_symbols; i++) {
      real_part = (bits(2*i)==0) ? M_1_DIV_SQRT2 : -M_1_DIV_SQRT2;
      imag_part = (bits(2*i+1)==0) ? M_1_DIV_SQRT2 : -M_1_DIV_SQRT2;
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
    int i, no_symbols = signal.size();
    out.set_size(2*no_symbols, false);
    for (i=0; i<no_symbols; i++) {
      out(2*i)   = ((std::real(signal(i)) > 0) ? 0 : 1);
      out(2*i+1) = ((std::imag(signal(i)) > 0) ? 0 : 1);
    }
  }

  bvec QPSK::demodulate_bits(const cvec &signal) const
  {
    bvec out;
    demodulate_bits(signal, out);
    return out;
  }


  //Outputs log-likelihood ratio of log (Pr(bit = 0|rx_symbols)/Pr(bit = 1|rx_symbols))
  void QPSK::demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits) const
  {
    soft_bits.set_size(2*rx_symbols.size(), false);
    double factor = 2*sqrt(2.0)/N0;
    for (int i=0; i<rx_symbols.size(); i++) {
      soft_bits(2*i) = std::real(rx_symbols(i))*factor;
      soft_bits(2*i+1) = std::imag(rx_symbols(i))*factor;
    }
  }

  //Outputs log-likelihood ratio for fading channels
  void QPSK::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const
  {
    soft_bits.set_size(2*rx_symbols.size(), false);
    std::complex<double> temp;
    double factor = 2*sqrt(2.0)/N0;
    
    for (int i=0; i<rx_symbols.size(); i++) {
      temp = rx_symbols(i)*std::conj(channel(i));
      soft_bits(2*i) = std::real(temp)*factor;
      soft_bits(2*i+1) = std::imag(temp)*factor;
    }
  }


  void QPSK::demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits) const
  {
    demodulate_soft_bits(rx_symbols, N0, soft_bits);
  }

  void QPSK::demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const
  {
    demodulate_soft_bits(rx_symbols, channel, N0, soft_bits);
  }


  //------------- class: PSK ----------------

  void PSK::modulate_bits(const bvec &bits, cvec &out) const
  {
    int no_symbols, i, symb;
    int bits_len = bits.length();

    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits_len % k) {
      it_warning("PSK::modulate_bits(): The number of input bits is not a multiple of log2(M),\nwhere M is a constellation size. Remainder bits are not modulated.");
    }
    no_symbols = floor_i(double(bits_len) / double(k));

    out.set_size(no_symbols, false);

    for (i=0; i<no_symbols; i++) {
      symb = bin2dec(bits.mid(i*k,k));
      out(i) = symbols(bits2symbols(symb));
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
    int i, est_symbol;
    out.set_size(k*signal.size(), false);
    double ang, temp;

    for (i=0; i<signal.size(); i++) {
      ang = std::arg(signal(i));
      temp = ( (ang < 0) ? (2*pi+ang) : ang );
      est_symbol = round_i(temp*(M/2)/pi) % M;
      out.replace_mid(i*k, bitmap.get_row(est_symbol));
    }
  }

  bvec PSK::demodulate_bits(const cvec &signal) const
  {
    bvec temp;
    demodulate_bits(signal, temp);
    return temp;
  }

  void PSK::demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double P0, P1;
    double d0min, d0min2, d1min, d1min2;
    double treshhold = -log(eps)*N0; // To be sure that any precision is left in the calculatation
    vec d(M), expd(M);
    double inv_N0 = 1/N0;

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {

      for (j=0; j<M; j++) {
	d(j) = sqr(rx_symbols(l)-symbols(j));
	expd(j) = exp((-d(j))*inv_N0);
      }

      for (i=0; i<k; i++) {
	P0 = P1 = 0;
	d0min = d0min2 = d1min = d1min2 = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min2 = d0min; d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min2 = d1min; d1min = d(S1(i,j)); }
	  P0 += expd(S0(i,j));
	  P1 += expd(S1(i,j));  
	}
	if ( (d0min2-d0min) > treshhold && (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
	else if ( (d0min2-d0min) > treshhold )
	  soft_bits(l*k+i) = -d0min*inv_N0-log(P1);
	else if ( (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = log(P0) + d1min*inv_N0;
	else
	  soft_bits(l*k+i) = log(P0/P1);
      }
    }
  }

  void PSK::demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double d0min, d1min;
    double inv_N0 = 1/N0;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {

      for (j=0; j<M; j++)
	d(j) = sqr(rx_symbols(l)-symbols(j));

      for (i=0; i<k; i++) {
	d0min = d1min = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min = d(S1(i,j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
      }
    }
 }

  void PSK::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double P0, P1;
    double d0min, d0min2, d1min, d1min2;
    double treshhold = -log(eps)*N0; // To be sure that any precision is left in the calculatation
    double inv_N0 = 1/N0;
    vec d(M), expd(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {

      for (j=0; j<M; j++) {
	d(j) = sqr(rx_symbols(l)-channel(l)*symbols(j));
	expd(j) = exp((-d(j))*inv_N0);
      }

      for (i=0; i<k; i++) {
	P0 = P1 = 0;
	d0min = d0min2 = d1min = d1min2 = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min2 = d0min; d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min2 = d1min; d1min = d(S1(i,j)); }
	  P0 += expd(S0(i,j));
	  P1 += expd(S1(i,j)); 
	}
	if ( (d0min2-d0min) > treshhold && (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
	else if ( (d0min2-d0min) > treshhold )
	  soft_bits(l*k+i) = -d0min*inv_N0-log(P1);
	else if ( (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = log(P0) + d1min*inv_N0;
	else
	  soft_bits(l*k+i) = log(P0/P1);
      }
    }
  }

  void PSK::demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double d0min, d1min;
    double inv_N0 = 1/N0;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {

      for (j=0; j<M; j++)
	d(j) = sqr(rx_symbols(l)-channel(l)*symbols(j));

      for (i=0; i<k; i++) {
	d0min = d1min = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min = d(S1(i,j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
      }
    }
  }

  void PSK::set_M(int Mary)
  {
    k = round_i(log2(double(Mary)));
    M = Mary;
    it_error_if(abs(pow2i(k)-Mary)>0.0,"M-ary PSK: M is not a power of 2");
    symbols.set_size(M, false);
    bits2symbols.set_size(M, false);
    bitmap = graycode(k);

    int i;
    double delta = 2.0*pi/M, epsilon = delta/10000.0;
    std::complex<double> symb;
    for (i=0; i<M; i++) {
      symb = std::complex<double>(std::polar(1.0,delta*i));
      if (std::abs(std::real(symb)) < epsilon) { symbols(i) = std::complex<double>(0.0,std::imag(symb)); }
      else if (std::abs(std::imag(symb)) < epsilon) { symbols(i) = std::complex<double>(std::real(symb),0.0); }
      else { symbols(i) = symb; }

      bits2symbols(bin2dec(bitmap.get_row(i))) = i;
    }

    //Calculate the soft bit mapping matrices S0 and S1
    S0.set_size(k,M/2,false);
    S1.set_size(k,M/2,false);
    int count0, count1, kk, m;
    bvec bits;
  
    for (kk=0; kk<k; kk++) {
      count0 = 0; 
      count1 = 0;
      for (m=0; m<M; m++) {
	bits = bitmap.get_row(m);
	if (bits(kk)==bin(0)) {
	  S0(kk,count0) = m;
	  count0++; 
	} else {
	  S1(kk,count1) = m;
	  count1++;
	}
      }
    }
  }

  //------------- class: QAM ----------------

  void QAM::modulate_bits(const bvec &bits, cvec &out) const
  {
    int no_symbols, i, symb;
    int bits_len = bits.length();

    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits_len % k) {
      it_warning("QAM::modulate_bits(): The number of input bits is not a multiple of log2(M),\nwhere M is a constellation size. Remainder bits are not modulated.");
    }
    no_symbols = floor_i(double(bits_len) / double(k));

    out.set_size(no_symbols, false);

    for (i=0; i<no_symbols; i++) {
      symb = bin2dec(bits.mid(i*k,k));
      out(i) = symbols(bits2symbols(symb));
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
    int i;
    out.set_size(k*signal.size(), false);

    int temp_real, temp_imag;

    for (i=0; i<signal.size(); i++) {
      temp_real = round_i((L-1)-(std::real(signal(i)*scaling_factor)+(L-1))/2.0);
      temp_imag = round_i((L-1)-(std::imag(signal(i)*scaling_factor)+(L-1))/2.0);
      if (temp_real<0) temp_real = 0; else if (temp_real>(L-1)) temp_real = (L-1);
      if (temp_imag<0) temp_imag = 0; else if (temp_imag>(L-1)) temp_imag = (L-1);
      out.replace_mid( k*i, bitmap.get_row(temp_imag*L + temp_real) );
    }
  }

  bvec QAM::demodulate_bits(const cvec &signal) const
  {
    bvec temp;
    demodulate_bits(signal, temp);
    return temp;
  }

  void QAM::demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double P0, P1;
    double d0min, d0min2, d1min, d1min2;
    double inv_N0 = 1/N0;
    double treshhold = -log(eps)*N0; // To be sure that any precision is left in the calculatation
    vec d(M), expd(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {

      for (j=0; j<M; j++) {
	d(j) = sqr(rx_symbols(l)-symbols(j));
	expd(j) = exp((-d(j))*inv_N0);
      }

      for (i=0; i<k; i++) {
	P0 = P1 = 0;
	d0min = d0min2 = d1min = d1min2 = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min2 = d0min; d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min2 = d1min; d1min = d(S1(i,j)); }
	  P0 += expd(S0(i,j));
	  P1 += expd(S1(i,j));  
	}
	if ( (d0min2-d0min) > treshhold && (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
	else if ( (d0min2-d0min) > treshhold )
	  soft_bits(l*k+i) = -d0min*inv_N0-log(P1);
	else if ( (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = log(P0) + d1min*inv_N0;
	else
	  soft_bits(l*k+i) = log(P0/P1);
      }
    }
  }


  void QAM::demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double d0min, d1min;
    double inv_N0 = 1/N0;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {

      for (j=0; j<M; j++)
	d(j) = sqr(rx_symbols(l)-symbols(j));

      for (i=0; i<k; i++) {
	d0min = d1min = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min = d(S1(i,j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
      }
    }
  }

  void QAM::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double P0, P1;
    double d0min, d0min2, d1min, d1min2;
    double treshhold = -log(eps)*N0; // To be sure that any precision is left in the calculatation
    double inv_N0 = 1/N0;
    vec d(M), expd(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {

      for (j=0; j<M; j++) {
	d(j) = sqr(rx_symbols(l)-channel(l)*symbols(j));
	expd(j) = exp((-d(j))*inv_N0);
      }

      for (i=0; i<k; i++) {
	P0 = P1 = 0;
	d0min = d0min2 = d1min = d1min2 = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min2 = d0min; d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min2 = d1min; d1min = d(S1(i,j)); }
	  P0 += expd(S0(i,j));
	  P1 += expd(S1(i,j)); 
	}
	if ( (d0min2-d0min) > treshhold && (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
	else if ( (d0min2-d0min) > treshhold )
	  soft_bits(l*k+i) = -d0min*inv_N0-log(P1);
	else if ( (d1min2-d1min) > treshhold )
	  soft_bits(l*k+i) = log(P0) + d1min*inv_N0;
	else
	  soft_bits(l*k+i) = log(P0/P1);
      }
    }
  }


  void QAM::demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const
  {
    int l, i, j;
    double d0min, d1min;
    double inv_N0 = 1/N0;
    vec d(M);

    soft_bits.set_size(k*rx_symbols.size(), false);

    for (l=0; l<rx_symbols.size(); l++) {

      for (j=0; j<M; j++)
	d(j) = sqr(rx_symbols(l)-channel(l)*symbols(j));

      for (i=0; i<k; i++) {
	d0min = d1min = 1e100;
	for (j=0; j<(M/2); j++) {
	  if (d(S0(i,j)) < d0min) { d0min = d(S0(i,j)); }
	  if (d(S1(i,j)) < d1min) { d1min = d(S1(i,j)); }
	}
	soft_bits(l*k+i) = (-d0min + d1min)*inv_N0;
      }
    }
  }


  void QAM::set_M(int Mary)
  {
    k = round_i(log2(double(Mary)));
    M = Mary;
    L = round_i(sqrt((double)M));
    it_error_if(abs(pow2i(k)-Mary)>0.01,"M-ary QAM: M is not a power of 2");

    int i, j, kk, m, count0, count1;
    bvec bits;

    symbols.set_size(M, false);
    bitmap.set_size(M, k, false);
    bits2symbols.set_size(M, false);
    bmat gray_code=graycode(round_i(log2(double(L))));
    average_energy = double(M-1)*2.0/3.0;
    scaling_factor = sqrt(average_energy);

    for (i=0; i<L; i++) {
      for (j=0; j<L; j++) {
	symbols(i*L+j) = std::complex<double>( ((L-1)-j*2)/scaling_factor ,((L-1)-i*2)/scaling_factor);
	bitmap.set_row( i*L+j, concat(gray_code.get_row(i), gray_code.get_row(j)) );
	bits2symbols( bin2dec(bitmap.get_row(i*L+j)) ) = i*L+j;
      }
    }

    //Calculate the soft bit mapping matrices S0 and S1
    S0.set_size(k,M/2,false);
    S1.set_size(k,M/2,false);
  
    for (kk=0; kk<k; kk++) {
      count0 = 0; 
      count1 = 0;
      for (m=0; m<M; m++) {
	bits = bitmap.get_row(m);
	if (bits(kk)==bin(0)) {
	  S0(kk,count0) = m;
	  count0++; 
	} else {
	  S1(kk,count1) = m;
	  count1++;
	}
      }
    }
  }

} //namespace itpp
