/*!
 * \file 
 * \brief One- and two-dimensional modulators - source file
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
#include <itpp/base/specmat.h>
#include <itpp/base/converters.h>
#include <itpp/base/logexpfunc.h>
#include <itpp/base/elmatfunc.h>


namespace itpp {


  // ----------------------------------------------------------------------
  // Modulator_1D (real domain signals)
  // ----------------------------------------------------------------------

  Modulator_1D::Modulator_1D():
    setup_done(false), k(0), M(0), bitmap(""), bits2symbols(""), symbols(""), 
    S0(""), S1("") {}

  Modulator_1D::Modulator_1D(const vec &symbols, const ivec &bits2symbols)
  {
    set(symbols, bits2symbols);
  }


  void Modulator_1D::set(const vec &in_symbols, const ivec &in_bits2symbols)
  {
    it_assert(in_symbols.size() == in_bits2symbols.size(), "Modulator_1D::set(): Number of symbols and bitmap does not match");
    M = in_bits2symbols.size();
    k = needed_bits(M - 1);
    symbols = in_symbols;
    bits2symbols = in_bits2symbols; 
    bitmap.set_size(M, k);
    for (int m = 0; m < bits2symbols.size(); m++) {
      bitmap.set_row(m, dec2bin(k, bits2symbols(m)));
    }
    calculate_softbit_matricies(bits2symbols);
    setup_done = true;
  }


  vec Modulator_1D::modulate(const ivec &symbolnumbers) const
  {
    it_assert0(setup_done, "Modulator_1D::modulate(): Modulator not ready.");
    vec temp(symbolnumbers.size());
    for (int i = 0; i < symbolnumbers.size(); i++)
      temp(i) = symbols(symbolnumbers(i));
    return temp;
  }

  ivec Modulator_1D::demodulate(const vec &signal) const
  {
    it_assert0(setup_done, "Modulator_1D::demodulate(): Modulator not ready.");
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


  void Modulator_1D::modulate_bits(const bvec &bits, vec& out) const
  {
    it_assert0(setup_done, "Modulator_1D::modulate_bits(): Modulator not ready.");
    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits.length() % k) {
      it_warning("Modulator_1D::modulate_bits(): The number of input bits is not a multiple of k (number of bits per symbol). Remainder bits are not modulated.");
    }
    int no_symbols = bits.length() / k;
    int pos;
    out.set_size(no_symbols);

    for (int i = 0; i < no_symbols; i++) {
      pos = 0;
      while (bits2symbols(pos) != bin2dec(bits.mid(i*k, k))) { pos++; }
      out(i) = symbols(pos);
    }
  }

  vec Modulator_1D::modulate_bits(const bvec &bits) const
  {
    vec out;
    modulate_bits(bits, out);
    return out;
  }


  void Modulator_1D::demodulate_bits(const vec &signal, bvec& out) const
  {
    it_assert0(setup_done, "Modulator_1D::demodulate_bits(): Modulator not ready.");
    double dist, mindist;
    int closest;
    out.set_size(k*signal.size());

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
      out.replace_mid(i*k, bitmap.get_row(closest));
    }
  }

  bvec Modulator_1D::demodulate_bits(const vec &signal) const
  {
    bvec out;
    demodulate_bits(signal, out);
    return out;
  }


  void Modulator_1D::calculate_softbit_matricies(const ivec& in_bits2symbols)
  {
    int count0, count1;

    // Allocate storage space for the result matricies:
    S0.set_size(k, M >> 1, false);
    S1.set_size(k, M >> 1, false);

    for (int i = 0; i < k; i++) {
      count0 = 0;
      count1 = 0;
      for (int j = 0; j < M; j++) {
	if (bitmap(j, i) == bin(0)) {
	  S0(i, count0++) = j;
	} 
	else {
	  S1(i, count1++) = j;
	}
      }
    }
  }



  // ----------------------------------------------------------------------
  // Modulator_2D (base class)
  // ----------------------------------------------------------------------

  Modulator_2D::Modulator_2D() : 
    setup_done(false), k(0), M(0), bitmap(""), bits2symbols(""), symbols(""),
    S0(""), S1("") {}

  Modulator_2D::Modulator_2D(const cvec &symbols, const ivec &bits2symbols)
  {
    set(symbols, bits2symbols);
  }

  void Modulator_2D::set(const cvec &in_symbols, const ivec &in_bits2symbols)
  {
    it_assert(in_symbols.size() == in_bits2symbols.size(), "Modulator_2D::set() Number of symbols and bits2symbols does not match");
    symbols = in_symbols;
    bits2symbols = in_bits2symbols; 
    M = bits2symbols.size();
    k = needed_bits(M - 1);
    bitmap.set_size(M, k);
    for (int m = 0; m < bits2symbols.size(); m++) {
      bitmap.set_row(m, dec2bin(k, bits2symbols(m)));
    }
    calculate_softbit_matricies(bits2symbols);
    setup_done = true;
  }


  void Modulator_2D::modulate(const ivec &symbolnumbers, cvec& out) const
  {
    it_assert0(setup_done, "Modulator_2D::modulate(): Modulator not ready.");
    out.set_size(symbolnumbers.length());
    for (int i = 0; i < symbolnumbers.length(); i++)
      out(i) = symbols(symbolnumbers(i));
  }

  cvec Modulator_2D::modulate(const ivec &symbolnumbers) const
  {
    cvec out(symbolnumbers.length());
    modulate(symbolnumbers, out);
    return out;
  }


  void Modulator_2D::demodulate(const cvec &signal, ivec& out) const
  {
    it_assert0(setup_done, "Modulator_2D::demodulate(): Modulator not ready.");
    double dist, mindist;
    int closest;
    out.set_size(signal.size());

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
      out(i) = closest;
    }
  }

  ivec Modulator_2D::demodulate(const cvec& signal) const
  {
    ivec out(signal.length());
    demodulate(signal, out);
    return out;
  }


  void Modulator_2D::modulate_bits(const bvec &bits, cvec &output) const
  {
    it_assert0(setup_done, "Modulator_2D::modulate_bits(): Modulator not ready.");
    // Check if some bits have to be cut and print warning message in such
    // case.
    if (bits.length() % k) {
      it_warning("Modulator_2D::modulate_bits(): The number of input bits is not a multiple of k (number of bits per symbol). Remainder bits are not modulated.");
    }
    int no_symbols = bits.length() / k;
    int pos;
    output.set_size(no_symbols);

    for (int i = 0; i < no_symbols; i++) {
      pos = 0;
      while (bits2symbols(pos) != bin2dec(bits.mid(i*k, k))) { pos++; }
      output(i) = symbols(pos);
    }
  }

  cvec Modulator_2D::modulate_bits(const bvec &bits) const
  {
    cvec output;
    modulate_bits(bits, output);
    return output;
  }

  void Modulator_2D::demodulate_bits(const cvec &signal, bvec &bits) const
  {
    it_assert0(setup_done, "Modulator_2D::demodulate_bist(): Modulator not ready.");
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
      bits.replace_mid(i*k, bitmap.get_row(closest));
    }
  }

  bvec Modulator_2D::demodulate_bits(const cvec &signal) const
  {
    bvec bits;
    demodulate_bits(signal, bits);
    return bits;
  }


  void Modulator_2D::demodulate_soft_bits(const cvec &rx_symbols, double N0, 
					  vec &soft_bits) const
  {
    it_assert0(setup_done, "Modulator_2D::demodulate_soft_bits(): Modulator not ready.");
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

  vec Modulator_2D::demodulate_soft_bits(const cvec &rx_symbols, 
					 double N0) const
  {
    vec out;
    demodulate_soft_bits(rx_symbols, N0, out);
    return out;
  }


  void Modulator_2D::demodulate_soft_bits(const cvec &rx_symbols, 
					  const cvec &channel, double N0, 
					  vec &soft_bits) const
  {
    it_assert0(setup_done, "Modulator_2D::demodulate_soft_bits(): Modulator not ready.");
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

  vec Modulator_2D::demodulate_soft_bits(const cvec &rx_symbols, 
					 const cvec &channel,
					 double N0) const
  {
    vec out;
    demodulate_soft_bits(rx_symbols, channel, N0, out);
    return out;
  }


  void Modulator_2D::demodulate_soft_bits_approx(const cvec &rx_symbols,
						 double N0, 
						 vec &soft_bits) const
  {
    it_assert0(setup_done, "Modulator_2D::demodulate_soft_bits_approx(): Modulator not ready.");
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

  vec Modulator_2D::demodulate_soft_bits_approx(const cvec &rx_symbols,
						double N0) const
  {
    vec out;
    demodulate_soft_bits_approx(rx_symbols, N0, out);
    return out;
  }


  void Modulator_2D::demodulate_soft_bits_approx(const cvec &rx_symbols, 
						 const cvec &channel, double N0,
						 vec &soft_bits) const
  {
    it_assert0(setup_done, "Modulator_2D::demodulate_soft_bits_approx(): Modulator not ready.");
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

  vec Modulator_2D::demodulate_soft_bits_approx(const cvec &rx_symbols, 
						const cvec &channel,
						double N0) const
  {
    vec out;
    demodulate_soft_bits_approx(rx_symbols, channel, N0, out);
    return out;
  }


  void Modulator_2D::calculate_softbit_matricies(const ivec& in_bits2symbols)
  {
    int count0, count1;

    // Allocate storage space for the result matricies:
    S0.set_size(k, M >> 1, false);
    S1.set_size(k, M >> 1, false);

    for (int i = 0; i < k; i++) {
      count0 = 0;
      count1 = 0;
      for (int j = 0; j < M; j++) {
	if (bitmap(j, i) == bin(0)) {
	  S0(i, count0++) = j;
	} 
	else {
	  S1(i, count1++) = j;
	}
      }
    }
  }


  // ----------------------------------------------------------------------
  // QAM
  // ----------------------------------------------------------------------

  void QAM::set_M(int Mary)
  {
    k = needed_bits(Mary - 1);
    M = Mary;
    it_assert(pow2i(k) == M, "QAM::set_M(): M is not a power of 2");
    L = round_i(std::sqrt(static_cast<double>(M)));

    double average_energy = (M - 1) * 2.0 / 3.0;
    scaling_factor = std::sqrt(average_energy);

    symbols.set_size(M);
    bitmap.set_size(M, k);
    bits2symbols.set_size(M);

    bmat gray_code = graycode(needed_bits(L - 1));

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
    calculate_softbit_matricies(bits2symbols);

    setup_done = true;
  }


  void QAM::modulate_bits(const bvec &bits, cvec &out) const
  {
    it_assert0(setup_done, "QAM::modulate_bits(): Modulator not ready.");
    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits.length() % k) {
      it_warning("QAM::modulate_bits(): The number of input bits is not a multiple of k (number of bits per symbol). Remainder bits are not modulated.");
    }
    int no_symbols = bits.length() / k;

    out.set_size(no_symbols, false);
    for (int i = 0; i < no_symbols; i++) {
      out(i) = symbols(bits2symbols(bin2dec(bits.mid(i*k, k))));
    }
  }

  cvec QAM::modulate_bits(const bvec &bits) const
  {
    cvec out;
    modulate_bits(bits, out);
    return out;
  }


  void QAM::demodulate_bits(const cvec &signal, bvec &out) const
  {
    it_assert0(setup_done, "QAM::demodulate_bits(): Modulator not ready.");
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
    bvec out;
    demodulate_bits(signal, out);
    return out;
  }


  // ----------------------------------------------------------------------
  // PSK
  // ----------------------------------------------------------------------

  void PSK::set_M(int Mary)
  {
    k = needed_bits(Mary - 1);
    M = Mary;
    it_assert(pow2i(k) == M, "PSK::set_M(): M is not a power of 2");

    symbols.set_size(M);
    bitmap = graycode(k);
    bits2symbols.set_size(M);

    double delta = m_2pi / M;
    double epsilon = delta / 10000.0;
    std::complex<double> symb;
    for (int i = 0; i < M; i++) {
      symb = std::complex<double>(std::polar(1.0, delta * i));
      if (std::fabs(std::real(symb)) < epsilon) { 
	symbols(i) = std::complex<double>(0.0, std::imag(symb)); 
      }
      else if (std::fabs(std::imag(symb)) < epsilon) { 
	symbols(i) = std::complex<double>(std::real(symb), 0.0);
      }
      else { 
	symbols(i) = symb; 
      }

      bits2symbols(bin2dec(bitmap.get_row(i))) = i;
    }

    // Calculate the soft bit mapping matrices S0 and S1
    calculate_softbit_matricies(bits2symbols);

    setup_done = true;
  }


  void PSK::modulate_bits(const bvec &bits, cvec &out) const
  {
    it_assert0(setup_done, "PSK::modulate_bits(): Modulator not ready.");
    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits.length() % k) {
      it_warning("PSK::modulate_bits(): The number of input bits is not a multiple of k (number of bits per symbol). Remainder bits are not modulated.");
    }
    int no_symbols = bits.length() / k;

    out.set_size(no_symbols, false);
    for (int i = 0; i < no_symbols; i++) {
      out(i) = symbols(bits2symbols(bin2dec(bits.mid(i*k, k))));
    }
  }

  cvec PSK::modulate_bits(const bvec &bits) const
  {
    cvec out;
    modulate_bits(bits, out);
    return out;
  }


  void PSK::demodulate_bits(const cvec &signal, bvec &out) const
  {
    it_assert0(setup_done, "PSK::demodulate_bits(): Modulator not ready.");
    int est_symbol;
    double ang, temp;

    out.set_size(k*signal.size(), false);

    for (int i = 0; i < signal.size(); i++) {
      ang = std::arg(signal(i));
      temp = (ang < 0) ? (m_2pi + ang) : ang;
      est_symbol = round_i(temp * (M >> 1) / pi) % M;
      out.replace_mid(i*k, bitmap.get_row(est_symbol));
    }
  }

  bvec PSK::demodulate_bits(const cvec &signal) const
  {
    bvec out;
    demodulate_bits(signal, out);
    return out;
  }


  // ----------------------------------------------------------------------
  // QPSK : PSK : Modulator_2D
  // ----------------------------------------------------------------------

  // These two specialised functions are defined for 4-QAM, not QPSK
  // modulator!!! 

//   void QPSK::demodulate_soft_bits(const cvec &rx_symbols, double N0,
// 				  vec &soft_bits) const
//   {
//     soft_bits.set_size(2*rx_symbols.size(), false);
//     double factor = 2 * std::sqrt(2.0) / N0;
//     for (int i = 0; i < rx_symbols.size(); i++) {
//       soft_bits(2*i) = std::real(rx_symbols(i)) * factor;
//       soft_bits(2*i+1) = std::imag(rx_symbols(i)) * factor;
//     }
//   }

//   void QPSK::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
// 				  double N0, vec &soft_bits) const
//   {
//     soft_bits.set_size(2*rx_symbols.size(), false);
//     std::complex<double> temp;
//     double factor = 2 * std::sqrt(2.0) / N0;
    
//     for (int i = 0; i < rx_symbols.size(); i++) {
//       temp = rx_symbols(i) * std::conj(channel(i));
//       soft_bits(2*i) = std::real(temp) * factor;
//       soft_bits(2*i+1) = std::imag(temp) * factor;
//     }
//   }


  // ----------------------------------------------------------------------
  // BPSK : PSK : Modulator_2D
  // ----------------------------------------------------------------------

  // Real signal:

  void BPSK::modulate_bits(const bvec &bits, vec &out) const
  {
    out.set_size(bits.size(), false);
    for (int i = 0; i < bits.size(); i++) { 
      out(i) = (bits(i) == 0 ? 1.0 : -1.0);
    }
  }

  void BPSK::demodulate_bits(const vec &signal, bvec &out) const
  {
    out.set_size(signal.size(), false);
    for (int i = 0; i < signal.length(); i++) { 
      out(i) = (signal(i) > 0) ? bin(0) : bin(1); 
    }
  }


  // Complex signal (imaginary part equal to 0):

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


  void BPSK::demodulate_bits(const cvec &signal, bvec &out) const
  {
    out.set_size(signal.size(), false);
    for (int i = 0; i < signal.length(); i++) { 
      out(i) = (std::real(signal(i)) > 0) ? bin(0) : bin(1); 
    }
  }

  bvec BPSK::demodulate_bits(const cvec &signal) const
  {
    bvec out(signal.size());
    demodulate_bits(signal, out);
    return out;
  }

  // Specialisation to decrease complexity

  void BPSK::demodulate_soft_bits(const cvec &rx_symbols, double N0,
				  vec &soft_bits) const
  {
    double factor = 4 / N0;
    soft_bits.set_size(rx_symbols.size(), false);

    for (int i = 0; i < rx_symbols.size(); i++) {
      soft_bits(i) = factor * std::real(rx_symbols(i));
    }
  }

  vec BPSK::demodulate_soft_bits(const cvec &rx_symbols, double N0) const
  {
    vec out;
    demodulate_soft_bits(rx_symbols, N0, out);
    return out;
  }


  void BPSK::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
				  double N0, vec &soft_bits) const
  {
    double factor = 4 / N0;
    soft_bits.set_size(rx_symbols.size(), false);

    for (int i = 0; i < rx_symbols.size(); i++) {
      soft_bits(i) = factor * std::real(rx_symbols(i) * std::conj(channel(i)));
    }
  }

  vec BPSK::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
				 double N0) const
  {
    vec out;
    demodulate_soft_bits(rx_symbols, channel, N0, out);
    return out;
  }


  void BPSK::demodulate_soft_bits_approx(const cvec &rx_symbols, double N0,
					 vec &soft_bits) const
  {
    demodulate_soft_bits(rx_symbols, N0, soft_bits);
  }

  vec BPSK::demodulate_soft_bits_approx(const cvec &rx_symbols, double N0) const
  {
    return demodulate_soft_bits(rx_symbols, N0);
  }


  void BPSK::demodulate_soft_bits_approx(const cvec &rx_symbols, 
					 const cvec &channel, double N0,
					 vec &soft_bits) const
  {
    demodulate_soft_bits(rx_symbols, channel, N0, soft_bits);
  }

  vec BPSK::demodulate_soft_bits_approx(const cvec &rx_symbols,
					const cvec &channel, double N0) const
  {
    return demodulate_soft_bits(rx_symbols, channel, N0);
  }


  // ----------------------------------------------------------------------
  // PAM : Modulator_2D
  // ----------------------------------------------------------------------

  void PAM::set_M(int Mary)
  {
    M = Mary;
    k = needed_bits(M - 1);
    it_assert(pow2i(k) == M, "PAM::set_M(): M is not a power of 2");

    symbols.set_size(M, false);
    bits2symbols.set_size(M, false);
    bitmap = graycode(k);
    double average_energy = (sqr(M) - 1) / 3.0;
    scaling_factor = std::sqrt(average_energy);

    for (int i = 0; i < M; i++) {
      symbols(i) = ((M-1) - i*2) / scaling_factor;
      bits2symbols(bin2dec(bitmap.get_row(i))) = i;
    }

    // Calculate the soft bit mapping matrices S0 and S1
    calculate_softbit_matricies(bits2symbols);

    setup_done = true;
  }


  // Real domain signal:

  void PAM::modulate_bits(const bvec &bits, vec &out) const
  {
    it_assert0(setup_done, "PAM::modulate_bits(): Modulator not ready.");
    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits.length() % k) {
      it_warning("PAM::modulate_bits(): The number of input bits is not a multiple of k (number of bits per symbol). Remainder bits are not modulated.");
    }
    int no_symbols = bits.length() / k;

    out.set_size(no_symbols, false);
    for (int i = 0; i < no_symbols; i++) {
      out(i) = std::real(symbols(bits2symbols(bin2dec(bits.mid(i*k, k)))));
    }
  }

  void PAM::demodulate_bits(const vec &signal, bvec &out) const
  {
    it_assert0(setup_done, "PAM::demodulate_bits(): Modulator not ready.");
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


  // Complex domain signal:

  void PAM::modulate_bits(const bvec &bits, cvec &out) const
  {
    it_assert0(setup_done, "PAM::modulate_bits(): Modulator not ready.");
    // Check if some bits have to be cut and print warning message in
    // such case.
    if (bits.length() % k) {
      it_warning("PAM::modulate_bits(): The number of input bits is not a multiple of k (number of bits per symbol). Remainder bits are not modulated.");
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

  void PAM::demodulate_bits(const cvec &signal, bvec &out) const
  {
    it_assert0(setup_done, "PAM::demodulate_bits(): Modulator not ready.");
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
    it_assert0(setup_done, "PAM::demodulate_soft_bits(): Modulator not ready.");
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

  vec PAM::demodulate_soft_bits(const cvec &rx_symbols, double N0) const
  {
    vec out;
    demodulate_soft_bits(rx_symbols, N0, out);
    return out;
  }


  void PAM::demodulate_soft_bits_approx(const cvec &rx_symbols, double N0,
					vec &soft_bits) const
  {
    it_assert0(setup_done, "PAM::demodulate_soft_bits_approx(): Modulator not ready.");
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

  vec PAM::demodulate_soft_bits_approx(const cvec &rx_symbols, double N0) const
  {
    vec out;
    demodulate_soft_bits_approx(rx_symbols, N0, out);
    return out;
  }


  void PAM::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
				 double N0, vec &soft_bits) const
  {
    it_assert0(setup_done, "PAM::demodulate_soft_bits(): Modulator not ready.");
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

  vec PAM::demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel,
				double N0) const
  {
    vec out;
    demodulate_soft_bits(rx_symbols, channel, N0, out);
    return out;
  }


  void PAM::demodulate_soft_bits_approx(const cvec &rx_symbols, 
					const cvec &channel, double N0,
					vec &soft_bits) const
  {
    it_assert0(setup_done, "PAM::demodulate_soft_bits_approx(): Modulator not ready.");
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

  vec PAM::demodulate_soft_bits_approx(const cvec &rx_symbols,
				       const cvec &channel,
				       double N0) const
  {
    vec out;
    demodulate_soft_bits_approx(rx_symbols, channel, N0, out);
    return out;
  }


} // namespace itpp
