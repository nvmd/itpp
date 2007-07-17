/*!
 * \file
 * \brief Implementation of a BCH encoder/decoder class
 * \author Pal Frenger, Steve Peters and Adam Piatyszek
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

#include <itpp/comm/bch.h>
#include <itpp/base/binary.h>
#include <itpp/base/specmat.h>

namespace itpp {

  //---------------------- BCH -----------------------------------

  BCH::BCH(int in_n, int in_k, int in_t, ivec genpolynom, bool sys) :
    n(in_n), k(in_k), t(in_t), systematic(sys)
  {
    //fix the generator polynomial g(x).
    ivec exponents(n-k+1);
    bvec temp = oct2bin(genpolynom);
    for (int i = 0; i < temp.length(); i++) {
      exponents(i) = static_cast<int>(temp(temp.length()-i-1)) - 1;
    }
    g.set(n+1, exponents);
  }

  void BCH::encode(const bvec &uncoded_bits, bvec &coded_bits)
  {
    int i, j, degree,
      itterations = floor_i(static_cast<double>(uncoded_bits.length()) / k);
    GFX m(n+1, k);
    GFX c(n+1, n);
    GFX r(n+1, n-k);
    GFX uncoded_shifted(n+1, n);
    coded_bits.set_size(itterations*n, false);
    bvec mbit(k), cbit(n);

    if (systematic)
      for (i = 0; i < n-k; i++)
	uncoded_shifted[i] = GF(n+1, -1);

    for (i = 0; i < itterations; i++) {
      //Fix the message polynom m(x).
      mbit = uncoded_bits.mid(i*k,k);
      for (j=0; j<k; j++) {
	degree = static_cast<int>(mbit(j)) - 1;
	m[j] = GF(n+1, degree);
	if (systematic) {
	  c[j] = m[j];
	  uncoded_shifted[j+n-k] = m[j];
	}
      }
      //Fix the outputbits cbit.
      if (systematic) {
	r = modgfx(uncoded_shifted, g);
	for (j = k; j < n; j++) {
	  c[j] = r[j-k];
	}
      } else {
	c = g * m;
      }
      for (j = 0; j < n; j++) {
	if (c[j] == GF(n+1, 0)) {
	  cbit(j) = 1;
	} else {
	  cbit(j) = 0;
	}
      }
      coded_bits.replace_mid(i*n, cbit);
    }
  }

  bvec BCH::encode(const bvec &uncoded_bits)
  {
    bvec coded_bits;
    encode(uncoded_bits, coded_bits);
    return coded_bits;
  }

  void BCH::decode(const bvec &coded_bits, bvec &decoded_bits)
  {
    int j, i, degree, kk, foundzeros, cisvalid,
      itterations = floor_i(static_cast<double>(coded_bits.length()) / n);
    bvec rbin(n), mbin(k);
    decoded_bits.set_size(itterations*k, false);

    GFX r(n+1, n-1), c(n+1, n-1), m(n+1, k-1), S(n+1, 2*t), Lambda(n+1),
      OldLambda(n+1), T(n+1), Ohmega(n+1), One(n+1, (char*)"0");
    GF delta(n+1), temp(n+1);
    ivec errorpos;

    for (i = 0; i < itterations; i++) {
      //Fix the received polynomial r(x)
      rbin = coded_bits.mid(i*n, n);
      for (j = 0; j < n; j++) {
	degree = static_cast<int>(rbin(j)) - 1;
	r[j] = GF(n+1, degree);
      }
      //Fix the syndrome polynomial S(x).
      S[0] = GF(n+1, -1);
      for (j = 1; j <= 2*t; j++) {
	S[j] =  r(GF(n+1, j));
      }
      if (S.get_true_degree() >= 1) { //Errors in the received word
	//Itterate to find Lambda(x).
	kk = 0;
	Lambda = GFX(n+1, (char*)"0");
	T = GFX(n+1, (char*)"0");
	while (kk < t) {
	  Ohmega = Lambda * (S + One);
	  delta = Ohmega[2*kk+1];
	  OldLambda = Lambda;
	  Lambda = OldLambda + delta*(GFX(n+1, (char*)"-1 0")*T);
	  if ((delta == GF(n+1,-1)) || (OldLambda.get_degree() > kk)) {
	    T = GFX(n+1, (char*)"-1 -1 0") * T;
	  } else {
	    T = (GFX(n+1, (char*)"-1 0") * OldLambda) / delta;
	  }
	  kk = kk + 1;
	}
	//Find the zeros to Lambda(x).
	errorpos.set_size(Lambda.get_true_degree(), true);
	foundzeros = 0;
	for (j = 0; j <= n-1; j++) {
	  temp = Lambda(GF(n+1, j));
	  if (Lambda(GF(n+1, j)) == GF(n+1, -1)) {
	    errorpos(foundzeros) = (n-j) % n;
	    foundzeros += 1;
	    if (foundzeros >= Lambda.get_true_degree()) {
	      break;
	    }
	  }
	}
	//Correct the codeword.
	for (j = 0; j < foundzeros; j++) {
	  rbin(errorpos(j)) += 1;
	}
	//Reconstruct the corrected codeword.
	for (j = 0; j < n; j++) {
	  degree = static_cast<int>(rbin(j)) - 1;
	  c[j] = GF(n+1, degree);
	}
	//Code word validation.
	S[0] = GF(n+1, -1);
	for (j = 1; j <= 2*t; j++) {
	  S[j] =  c(GF(n+1, j));
	}
	if (S.get_true_degree() <= 0) { //c(x) is a valid codeword.
	  cisvalid = true;
	} else {
	  cisvalid = false;
	}
      } else {
	c = r;
	cisvalid = true;
      }
      //Construct the message bit vector.
      if (cisvalid) { //c(x) is a valid codeword.
	if (c.get_true_degree() > 1) {
	  if (systematic) {
	    for (j = 0; j < k; j++)
	      m[j] = c[j];
	  } else {
	    m = divgfx(c, g);
	  }
	  mbin.clear();
	  for (j = 0; j <= m.get_true_degree(); j++) {
	    if (m[j] == GF(n+1, 0)) {
	      mbin(j) = 1;
	    }
	  }
	} else { //The zero word was transmitted
	  mbin = zeros_b(k);
	  m = GFX(n+1, (char*)"-1");
	}
      } else { //Decoder failure.
	mbin = zeros_b(k);
	m = GFX(n+1, (char*)"-1");
      }
      decoded_bits.replace_mid(i*k, mbin);
    }
  }


  bvec BCH::decode(const bvec &coded_bits)
  {
    bvec decoded_bits;
    decode(coded_bits, decoded_bits);
    return decoded_bits;
  }


  // --- Soft-decision decoding is not implemented ---

  void BCH::decode(const vec &received_signal, bvec &output)
  {
    it_error("BCH::decode(): Soft-decision decoding is not implemented");
  }

  bvec BCH::decode(const vec &received_signal)
  {
    it_error("BCH::decode(): Soft-decision decoding is not implemented");
    return bvec();
  }


} // namespace itpp
