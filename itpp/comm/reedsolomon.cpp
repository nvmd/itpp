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
  \brief Implementation of a Reed-Solomon codec class.
  \author Pål Frenger

  $Revision$

  $Date$
*/

#include <itpp/base/binary.h>
#include <itpp/comm/reedsolomon.h>

namespace itpp { 

  //-------------------- Help Function ----------------------------

  //! Local help function
  GFX formal_derivate(const GFX &f)  {
    int degree = f.get_true_degree();
    int q = f.get_size();
    int i;
    GFX fprim(q,degree);
    fprim.clear();
    for (i=0; i<=degree-1; i+=2) {
      fprim[i] = f[i+1];
    }
    return fprim;
  }

  //-------------------- Reed-Solomon ----------------------------
  //A Reed-Solomon code is a q^m-ary BCH code of length n = pow(q,m)-1. 
  //k = pow(q,m)-1-t. This class works for q==2.
  Reed_Solomon::Reed_Solomon(int in_m, int in_t) {
    m = in_m;
    t = in_t;
    n = round_i(pow(2,m)-1);
    k = round_i(pow(2,m))-1-2*t;
    q = round_i(pow(2,m));
    int i;
    GFX x(q,(char *)"-1 0");
    ivec alphapow(1);
    g.set(q,(char *)"0");
    for(i=1; i<=2*t; i++) {
      alphapow(0) = i;
      g *= (x-GFX(q,alphapow));
    }
  }

  void Reed_Solomon::encode(const bvec &uncoded_bits, bvec &coded_bits)
  {
    int i, j, itterations = (int)floor( (double)uncoded_bits.length() / (k*m) );
    GFX mx(q,k), cx(q,n);
    GF mpow;
    bvec mbit(k*m), cbit(m);

    coded_bits.set_size(itterations*n*m, false);

    for (i=0; i<itterations; i++) {
      //Fix the message polynom m(x).
      for (j=0; j<k; j++) {
	mpow.set(q,uncoded_bits.mid((i*m*k)+(j*m),m));
	mx[j] = mpow;
      }  
      //Fix the outputbits cbit.
      cx = g*mx;
      for (j=0; j<n; j++) {
	cbit = cx[j].get_vectorspace();
	coded_bits.replace_mid((i*n*m)+(j*m), cbit);
      }
    }
  }

  bvec Reed_Solomon::encode(const bvec &uncoded_bits)
  {
    bvec coded_bits;
    encode(uncoded_bits, coded_bits);
    return coded_bits;
  }

  void Reed_Solomon::decode(const bvec &coded_bits, bvec &decoded_bits)
  {
    int j, i, kk, l, L, foundzeros, itterations = (int)floor( (double)coded_bits.length() / (n*m) ), decoderfailure;
    bvec mbit(m*k);
    decoded_bits.set_size(itterations*k*m, false);

    GFX rx(q,n-1), cx(q,n-1), mx(q,k-1), ex(q,n-1), S(q,2*t), Lambda(q), Lambdaprim(q), OldLambda(q), T(q), Ohmega(q);
    GFX dummy(q), One(q,(char*)"0"), Ohmegatemp(q);
    GF delta(q), tempsum(q), rtemp(q), temp(q), Xk(q), Xkinv(q);
    ivec errorpos;

    for (i=0; i<itterations; i++) {
      decoderfailure = false;
      //Fix the received polynomial r(x)
      for (j=0; j<n; j++) {  rtemp.set(q,coded_bits.mid(i*n*m + j*m, m)); rx[j] = rtemp; }
      //Fix the syndrome polynomial S(x).
      S.clear();
      for (j=1; j<=2*t; j++) { S[j] =  rx(GF(q,j)); }
      if (S.get_true_degree() == 0) {cx = rx; decoderfailure = false; }
      else {//Errors in the received word
	    //Itterate to find Lambda(x).
	kk = 0; Lambda = GFX(q,(char*)"0"); L = 0; T = GFX(q,(char*)"-1 0");   
	while (kk<2*t) {
	  kk = kk + 1;
	  tempsum = GF(q,-1);
	  for (l=1; l<=L; l++) { tempsum += Lambda[l] * S[kk-l]; }
	  delta = S[kk] - tempsum;
	  if (delta != GF(q,-1)) {
	    OldLambda = Lambda;
	    Lambda -= delta*T;
	    if (2*L<kk) { L = kk - L; T = OldLambda/delta;}
	  } 
	  T = GFX(q,(char*)"-1 0") * T;
	}
	//Find the zeros to Lambda(x).
	errorpos.set_size(Lambda.get_true_degree(), false);
	errorpos.clear();
	foundzeros = 0;
	for (j=q-2; j>=0; j--) {
	  temp = Lambda( GF(q,j) );
	  if  (Lambda( GF(q,j) ) == GF(q,-1) ) {
	    errorpos( foundzeros ) = (n-j) % n;
	    foundzeros +=1;
	    if (foundzeros >= Lambda.get_true_degree()) { break; }
	  }
	}
	if (foundzeros != Lambda.get_true_degree()) { decoderfailure = false; } 
	else {
	  //Compute Ohmega(x) using the key equation for RS-decoding
	  Ohmega.set_degree(2*t);
	  Ohmegatemp = Lambda * (One +S);
	  for (j=0; j<=2*t; j++) { Ohmega[j] = Ohmegatemp[j]; }
	  Lambdaprim = formal_derivate(Lambda);
	  //Find the error polynomial
	  ex.clear();
	  for (j=0; j<foundzeros; j++) {
	    Xk = GF(q,errorpos(j)); Xkinv = GF(q,0) / Xk;
	    ex[errorpos(j)] = (Xk * Ohmega(Xkinv)) / Lambdaprim(Xkinv);
	  }
	  //Reconstruct the corrected codeword.
	  cx = rx + ex;
	  //Code word validation
	  S.clear(); for (j=1; j<=2*t; j++) { S[j] =  rx(GF(q,j)); }
	  if (S.get_true_degree() >= 1) { decoderfailure=false; }
	}
      }
      //Find the message polynomial
      mbit.clear();
      if (decoderfailure == false) {
	if (cx.get_true_degree() >= 1) {// A nonzero codeword was transmitted
	  mx = divgfx(cx,g);
	  for (j=0; j<=mx.get_true_degree(); j++) { mbit.replace_mid(j*m,mx[j].get_vectorspace()); }
	}
      } 
      decoded_bits.replace_mid(i*m*k,mbit);
    }
  }

  bvec Reed_Solomon::decode(const bvec &coded_bits)
  {
    bvec decoded_bits;
    decode(coded_bits, decoded_bits);
    return decoded_bits;
  }

  // --------------- Soft-decision decoding is not implemented --------------------------------
  void Reed_Solomon::decode(const vec &received_signal, bvec &output)
  {
    it_error("Reed_Solomon::decode(vec, bvec); soft-decision decoding is not implemented");
  }

  bvec Reed_Solomon::decode(const vec &received_signal)
  {
    it_error("Reed_Solomon::decode(vec, bvec); soft-decision decoding is not implemented");
    return bvec();
  }


} //namespace itpp
