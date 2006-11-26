/*!
 * \file 
 * \brief Implementation of vector ("MIMO") modulator classes
 * \author  Erik G. Larsson
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

#include <itpp/comm/modulator_nd.h>
#include <itpp/base/converters.h>
#include <itpp/comm/commfunc.h>
#include <itpp/base/elmatfunc.h>
#include <itpp/base/cholesky.h>
#include <itpp/base/inv.h>

namespace itpp {

  // ----------------------- MIMO GENERAL -----------------------------
  
  QLLRvec Modulator_ND::probabilities(QLLR l)
  {
    QLLRvec result(2);
    
    if (l<0) {          // this can be done more efficiently 
      result(1) = -llrcalc.jaclog(0,-l);
      result(0) = result(1) - l;
    } else {
      result(0) = -llrcalc.jaclog(0,l);
      result(1) = result(0) + l;
    }
    return result;
  }
  
  Vec<QLLRvec> Modulator_ND::probabilities(QLLRvec &l)
  {
    Vec<QLLRvec> result(length(l));   
    for (int i=0; i<length(l); i++) {
      result(i) = probabilities(l(i));
    }
    return result;
  }

  void Modulator_ND::update_LLR(Vec<QLLRvec> &logP_apriori, QLLRvec &p1, QLLRvec &p0, short s, QLLR scaled_norm, short j)
  {
    QLLR log_apriori_prob_const_point = 0;
    short b=0;
    for (short i=0; i<k(j); i++) {
      log_apriori_prob_const_point += ((bitmap(j)(s,i)==0) ? logP_apriori(b)(1) : logP_apriori(b)(0));
      b++;
    }
    
    b=0;
    for (short i=0; i<k(j); i++) {
      if (bitmap(j)(s,i)==0) {
	p1(b) =  llrcalc.jaclog(p1(b), scaled_norm + log_apriori_prob_const_point );
      }  else {
	p0(b) = llrcalc.jaclog(p0(b),  scaled_norm + log_apriori_prob_const_point );
      }
      b++;
    }
  }

  void Modulator_ND::update_LLR(Vec<QLLRvec> &logP_apriori, QLLRvec &p1, QLLRvec &p0, svec &s, QLLR scaled_norm)
  {
    QLLR log_apriori_prob_const_point = 0;
    short b=0;
    for (short j=0; j<nt; j++) {
      for (short i=0; i<k(j); i++) {
	log_apriori_prob_const_point += ((bitmap(j)(s[j],i)==0) ? logP_apriori(b)(1) : logP_apriori(b)(0));
	b++;
      }
    }
    
    b=0;
    for (short j=0; j<nt; j++) {
      for (short i=0; i<k(j); i++) {
	if (bitmap(j)(s[j],i)==0) {
	  p1(b) =  llrcalc.jaclog(p1(b), scaled_norm + log_apriori_prob_const_point );
	}  else {
	  p0(b) = llrcalc.jaclog(p0(b),  scaled_norm + log_apriori_prob_const_point );
	}
	b++;
      }
    }	    
  }

  void Modulator_NRD::update_norm(double &norm, short k, short sold, short snew, vec &ytH, mat &HtH, svec &s)
  {
   
    short m = length(s);
    double cdiff = symbols(k)[snew]-symbols(k)[sold];
    
    norm += sqr(cdiff)*HtH(k,k);
    cdiff = 2.0*cdiff;
    norm -= cdiff*ytH[k];
     for (short i=0; i<m; i++) {
      norm += cdiff*HtH(i,k)*symbols(k)[s[i]];
     }
  }

  void Modulator_NCD::update_norm(double &norm, short k, short sold, short snew, cvec &ytH, cmat &HtH, svec &s)
  {
    short m = length(s);
    std::complex<double> cdiff = symbols(k)[snew]-symbols(k)[sold];
    
    norm += sqr(cdiff)*(HtH(k,k).real());
    cdiff = 2.0*cdiff;
    norm -= (cdiff.real()*ytH[k].real() - cdiff.imag()*ytH[k].imag());
    for (short i=0; i<m; i++) {
      norm += (cdiff*HtH(i,k)*conj(symbols(k)[s[i]])).real();
    }
  }

  void Modulator_NRD::map_demod(QLLRvec &LLR_apriori,  QLLRvec &LLR_aposteriori,  
			   double sigma2,  vec &h, vec &y)
  {
    it_assert(length(LLR_apriori)==sum(k),"Modulator_NRD::map_demod()");
    it_assert(length(LLR_apriori)==length(LLR_aposteriori),"Modulator_NRD::map_demod()");
    it_assert(length(h)==length(y),"Modulator_NRD::map_demod()");
    it_assert(length(h)==nt,"Modulator_NRD:map_demod()");

    int b=0;
    double oo_2s2 = 1.0/(2.0*sigma2);
    double norm2;
    QLLRvec temp, bnum, bdenom;
    for (short i=0; i<nt; i++) {
      temp=LLR_apriori(b,b+k(i)-1);
      bnum = (-QLLR_MAX)*ones_i(k(i));
      bdenom = (-QLLR_MAX)*ones_i(k(i));    
      Vec<QLLRvec> logP_apriori = probabilities(temp);
      for (short j=0; j<M(i); j++) {
	norm2 = sqr(y(i)-h(i)*symbols(i)(j));
	QLLR scaled_norm = llrcalc.to_qllr(-norm2*oo_2s2);
	update_LLR(logP_apriori, bnum, bdenom, j, scaled_norm,i);
      }
      LLR_aposteriori.set_subvector(b,bnum-bdenom);
      b+=k(i);
    }
  };

  void Modulator_NCD::map_demod(QLLRvec &LLR_apriori,  QLLRvec &LLR_aposteriori,  
				double sigma2,  cvec &h, cvec &y)
  {
    it_assert(length(LLR_apriori)==sum(k),"Modulator_NCD::map_demod()");
    it_assert(length(LLR_apriori)==length(LLR_aposteriori),"Modulator_NCD::map_demod()");
    it_assert(length(h)==length(y),"Modulator_NCD::map_demod()");
    it_assert(length(h)==nt,"Modulator_NCD:map_demod()");

    int b=0;
    double oo_s2 = 1.0/sigma2;
    double norm2;
    QLLRvec temp, bnum, bdenom;
    for (short i=0; i<nt; i++) {
      temp=LLR_apriori(b,b+k(i)-1);
      bnum = (-QLLR_MAX)*ones_i(k(i));
      bdenom = (-QLLR_MAX)*ones_i(k(i));    
      Vec<QLLRvec> logP_apriori = probabilities(temp);
      for (short j=0; j<M(i); j++) {
	norm2 = sqr(y(i)-h(i)*symbols(i)(j)); 
	QLLR scaled_norm = llrcalc.to_qllr(-norm2*oo_s2);
	update_LLR(logP_apriori, bnum, bdenom, j, scaled_norm,i);
      }
      LLR_aposteriori.set_subvector(b,bnum-bdenom);
      b+=k(i);
    }
  };

  void Modulator_NRD::map_demod(QLLRvec &LLR_apriori,  QLLRvec &LLR_aposteriori,  
			   double sigma2,  mat &H, vec &y)
  {
    short nr = H.rows();
    short np=sum(k); // number of bits in total
    it_assert(length(LLR_apriori)==np,"Modulator_NRD::map_demod()");
    it_assert(length(LLR_apriori)==length(LLR_aposteriori),"Modulator_NRD::map_demod()");
    it_assert(H.rows()==length(y),"Modulator_NRD::map_demod()");
    it_assert(H.cols()==nt,"Modulator_NRD:map_demod()");
    
    short mode=0;
     for (short i=0; i<length(M); i++) {
       if (nt*M(i)>4) { mode = 1; }    // differential update only pays off for larger dimensions
     }
    
    Vec<QLLRvec> logP_apriori = probabilities(LLR_apriori);
    
    mat HtH = H.transpose()*H; 
    vec ytH = H.transpose()*y; 
    
    QLLRvec bnum = (-QLLR_MAX)*ones_i(np);
    QLLRvec bdenom = (-QLLR_MAX)*ones_i(np);    
    svec s(nt);
    s.clear();
    double norm = 0.0;
    double oo_2s2 = 1.0/(2.0*sigma2);
    
    // Go over all constellation points  (r=dimension, s=vector of symbols)
    short r = nt-1;
    while (1==1) {
      
      if (mode==1) {
	update_norm(norm, r, s[r], 0, ytH, HtH, s);
      }
      s[r] = 0;      
      while (1==1) {
	if (s[r] > M(r)-1)  {
	  if (r==nt-1) {
	    goto exit_point;
	  }
	  r++;
	}   else {
	  if (r==0) {
	    if (mode==0) {
	      norm = 0.0;
	      for (short p=0; p<nr; p++) {
		double d = y[p];
		for (short i=0; i<nt; i++) {
		  d -= H(p,i)*symbols(i)[s[i]];
		}
		norm += sqr(d);
	      }
	    }
	    QLLR scaled_norm = llrcalc.to_qllr(-norm*oo_2s2);
	    update_LLR(logP_apriori, bnum, bdenom, s, scaled_norm);
	  } else {
	    r--;
	    break;
	  }
	}
	if (mode==1) {
	  update_norm(norm, r, s[r], s[r]+1, ytH, HtH, s);
	}
	s[r]++;
      }
    }
    
  exit_point:
    LLR_aposteriori = bnum - bdenom;

  };

  void Modulator_NCD::map_demod(QLLRvec &LLR_apriori,  QLLRvec &LLR_aposteriori,  double sigma2,  
				cmat &H, cvec &y)
  {
    short nr = H.rows();
    short np=sum(k); // number of bits in total
    it_assert(length(LLR_apriori)==np,"Modulator_NCD::map_demod()");
    it_assert(length(LLR_apriori)==length(LLR_aposteriori),"Modulator_NCD::map_demod()");
    it_assert(H.rows()==length(y),"Modulator_NCD::map_demod()");
    it_assert(H.cols()==nt,"Modulator_NCD:map_demod()");
    
    short mode=0;  
    for (short i=0; i<length(M); i++) {
      if (nt*M(i)>4) { mode = 1; }    // differential update only pays off for larger dimensions
    }
    
    Vec<QLLRvec> logP_apriori = probabilities(LLR_apriori);
    
    cmat HtH = H.hermitian_transpose()*H; 
    cvec ytH = conj(H.hermitian_transpose()*y); 
    
    QLLRvec bnum = (-QLLR_MAX)*ones_i(np);
    QLLRvec bdenom = (-QLLR_MAX)*ones_i(np);   
    svec s(nt);
    s.clear();
    double norm = 0.0;
    double oo_s2 = 1.0/sigma2;
    std::complex<double> d;

    // Go over all constellation points  (r=dimension, s=vector of symbols)
    short r = nt-1;
    while (1==1) {
      
      if (mode==1) {
	update_norm(norm, r, s[r], 0, ytH, HtH, s);
      }      
      s[r] = 0;
      while (1==1) {
	if (s[r] > M(r)-1)  {
	  if (r==nt-1) {
	    goto exit_point;
	  }
	  r++;
	}   else {
	  if (r==0) {
	    if (mode==0) {
	      norm = 0.0;
	      for (short p=0; p<nr; p++) {
		d = y[p];
		for (short i=0; i<nt; i++) {
		  d -= H(p,i)*symbols(i)[s[i]];
		}
		norm += sqr(d);
	      }
	    }
    	    QLLR scaled_norm = llrcalc.to_qllr(-norm*oo_s2);
	    update_LLR(logP_apriori, bnum, bdenom, s, scaled_norm);
	  } else {
	    r--;
	    break;
	  }
	}
	if (mode==1) {
	  update_norm(norm, r, s[r], s[r]+1, ytH, HtH, s);
	}
	s[r]++;
      }
    }
    
  exit_point:
    LLR_aposteriori = bnum - bdenom;
  };



  vec Modulator_NRD::modulate_bits(const bvec &bits) const
  {
    vec result(nt);

    it_assert(length(bits)==sum(k),"Modulator_NRD::modulate_bits(): The number of input bits does not match.");

    int b=0;
    for (short i=0; i<nt; i++) {
      short symb = bin2dec(bits.mid(b,k(i)));
      result(i) = symbols(i)(bits2symbols(i)(symb));
      b+=k(i);
    }

    return result;
  };


  cvec Modulator_NCD::modulate_bits(const bvec &bits) const
  {
    cvec result(nt);

    it_assert(length(bits)==sum(k),"Modulator_NCD::modulate_bits(): The number of input bits does not match.");

    int b=0;
    for (short i=0; i<nt; i++) {
      short symb = bin2dec(bits.mid(b,k(i)));
      result(i) = symbols(i)(bits2symbols(i)(symb));
      b+=k(i);
    }

    return result;
  };
  



  std::ostream &operator<<(std::ostream &os, const Modulator_NRD &mod)
  {
    os << "--- REAL MIMO (NRD) CHANNEL ---------" << std::endl;
    os << "Dimension (nt):           " << mod.nt << std::endl;
    os << "Bits per dimension (k):   " << mod.k << std::endl;
    os << "Symbols per dimension (M):" << mod.M << std::endl;
    for (int i=0; i<mod.nt; i++) { 
      os << "Bitmap for dimension " << i << ": " << mod.bitmap(i) << std::endl;
      // skip printing the trailing zero
      os << "Symbol coordinates for dimension " << i << ": " << mod.symbols(i).left(mod.M(i)) << std::endl;
    }
    os << mod.get_llrcalc() << std::endl;
    return os;
  };
  
  std::ostream &operator<<(std::ostream &os, const Modulator_NCD &mod)
  {
    os << "--- COMPLEX MIMO (NCD) CHANNEL --------" << std::endl;
    os << "Dimension (nt):           " << mod.nt << std::endl;
    os << "Bits per dimension (k):   " << mod.k << std::endl;
    os << "Symbols per dimension (M):" << mod.M << std::endl;
    for (int i=0; i<mod.nt; i++) { 
      os << "Bitmap for dimension " << i << ": " << mod.bitmap(i) << std::endl;
      os << "Symbol coordinates for dimension " << i << ": " << mod.symbols(i).left(mod.M(i)) << std::endl;
    }
    os << mod.get_llrcalc() << std::endl;
    return os;
  };

  // ------------------------- MIMO with uniform PAM ----------------------------

  ND_UPAM::ND_UPAM(int nt_in, int Mary)
  {
    nt=nt_in;
    set_Gray_PAM(nt,Mary);
  };
  
  void ND_UPAM::set_Gray_PAM(int nt_in, int Mary)
  {
    nt=nt_in;
    ivec Mary_temp(nt);
    for (int i=0; i<nt; i++) {
      Mary_temp(i)=Mary;
    }
    set_Gray_PAM(nt,Mary_temp);
  };

  void ND_UPAM::set_Gray_PAM(int nt_in, ivec Mary)
  {
    nt=nt_in;
    it_assert(length(Mary)==nt,"ND_UPAM::set_Gray_PAM() Mary has wrong length");
    k.set_size(nt);
    M=to_svec(Mary);
    bitmap.set_size(nt);
    symbols.set_size(nt);
    bits2symbols.set_size(nt);    
    spacing.set_size(nt);

    for (int i=0; i<nt; i++) {
      k(i) = round_i(log2(double(M(i))));
      it_assert( ((k(i)>0) && ((1<<k(i))==M(i))),"ND_UPAM::set_Gray_PAM(): M is not a power of 2.");
       
      symbols(i).set_size(M(i)+1);
      bits2symbols(i).set_size(M(i));
      bitmap(i) = graycode(k(i));
      double average_energy = double(M(i)*M(i)-1.0)/3.0;
      double scaling_factor = std::sqrt(average_energy);
      
      for (int j=0; j<M(i); j++) {
	symbols(i)(j) = double((M(i)-1)-j*2) / scaling_factor;
	bits2symbols(i)(bin2dec(bitmap(i).get_row(j))) = j;
      }      
      
      symbols(i)(M(i))=0.0;  // must end with a zero; only for a trick exploited in update_norm()

      spacing(i)=2.0/scaling_factor;
    }
  };
  

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#define sign(a)  (((a)>0)?(1):(-1))

  int ND_UPAM::sphere_search_SE(vec &y_in, mat &H, imat &zrange, double r, svec &zhat)
  {
    // The implementation of this function basically follows the
    // Schnorr-Eucner algorithm described in Agrell et al. (IEEE
    // Trans.  IT, 2002), but taking into account constellation
    // boundaries, see the "accelerated sphere decoder" in Boutros et
    // al. (IEEE Globecom, 2003).  No lattice reduction is performed.
    // Potentially the function can be speeded up by performing
    // lattice reduction, but it seems difficult to keep track of
    // constellation boundaries.
        
    mat R=chol(H.transpose()*H);
    mat Ri=inv(R);
    mat Q=H*Ri;
    vec y=Q.transpose()*y_in;
    mat Vi=Ri.transpose();

    short n=H.cols();
    vec dist(n);
    dist[n-1] = 0;
    double bestdist = r*r;
    short status = -1; // search failed 
    
    mat E=zeros(n,n);
    for (short i=0; i<n; i++) {   // E(k,:) = y*Vi; 
      for (short j=0; j<n; j++) {
	E(i*n+n-1) += y(j)*Vi(j+n*i);
      }
    }
    
    svec z(n);
    zhat.set_size(n);
    z(n-1) = (short) std::floor(0.5 + E(n*n-1));
    z(n-1) = max(z(n-1),zrange(n-1,0));
    z(n-1) = min(z(n-1),zrange(n-1,1));
    double p = (E(n*n-1)-z(n-1))/Vi(n*n-1);
    svec step(n);
    step(n-1) = sign(p);
    
    // Run search loop 
    short k=n-1;  // k uses natural indexing, goes from 0 to n-1 
    
    while (1==1) {
      double newdist = dist(k) + p*p;
      
      if ((newdist<bestdist) && (k!=0)) {
	
	for (short i=0; i<k; i++) {
	  E(k-1+i*n) = E(k+i*n) - p*Vi(k+i*n);
	}
	
	k--;
	dist(k) = newdist;
	z(k) = (short) std::floor(0.5+E(k+k*n));
	z(k) = max(z(k),zrange(k,0));
	z(k) = min(z(k),zrange(k,1));
	p = (E(k+k*n)-z(k))/Vi(k+k*n);
	
	step[k] = sign(p);
      } else {
	while (1==1) {
	  if (newdist<bestdist) {
	    for (short j=0; j<n; j++) { zhat(j) = z(j); }
	    bestdist = newdist; 
	    status = 0;
	  }
	  else if (k==n-1) {
	    goto exit_point; 
	  } else {
	    k++;
	  }
	  
	  z[k] += step(k);
	  
	  if ((z(k)<zrange(k,0)) || (z(k)>zrange(k,1))) {
	    step(k) = (-step(k) - sign(step(k)));
	    z(k) += step(k);
	  }
	  
	  if ((z(k)>=zrange(k,0)) && (z(k)<=zrange(k,1))) {
	    break;
	  }
	}
	
	p = (E(k+k*n)-z(k))/Vi(k+k*n);
	step(k) = (-step(k) - sign(step(k)));
      }
    }
    
  exit_point:
    return status;
  }
  
  
  int ND_UPAM::sphere_decoding(vec &y, mat &H, double rstart, double rmax, 
			       double stepup, QLLRvec &detected_bits)
  {
    it_assert(H.rows()==length(y),"ND_UPAM::sphere_decoding(): dimension mismatch");
    it_assert(H.rows()==nt,"ND_UPAM::sphere_decoding(): dimension mismatch");
    it_assert(rstart>0,"ND_UPAM::sphere_decoding(): radius error");
    it_assert(rmax>rstart,"ND_UPAM::sphere_decoding(): radius error");
    
    // This function can be improved, e.g., by using an ordered search.

    vec ytemp=y;
    mat Htemp(H.rows(),H.cols());
    for (short i=0; i<H.cols(); i++) {
      Htemp.set_col(i,H.get_col(i)*spacing(i));
      ytemp+=Htemp.get_col(i)*0.5*(M(i)-1.0);
    }

    imat crange(nt,2);
    for (short i=0; i<nt; i++) {
      crange(i,0) = 0;
      crange(i,1) = M(i)-1;
    }

    short status;  
    double r = rstart;
    svec s(sum(M));
    while (r<=rmax) {
      status = sphere_search_SE(ytemp,Htemp,crange,r,s);
      
      if (status==0) { // search successful
	detected_bits.set_size(sum(k));
	int b=0;
	for (short j=0; j<nt; j++) {
	  for (short i=0; i<k(j); i++) {
	    if (bitmap(j)((M(j)-1-s[j]),i)==0) {
	      detected_bits(b) = 1000;
	    }  else {
	      detected_bits(b) = -1000;
	    }
	    b++;
	  }
	}	    
	
	return status;
      }
      r = r*stepup;
    }

    return status;
    
  }
  

  // -------------------- MIMO with uniform QAM ---------------------------- 

  // The ND_UQAM class could alternatively have been implemented by
  // using a ND_UPAM class of twice the dimension, but this does not
  // fit as elegantly into the class structure

  ND_UQAM::ND_UQAM(int nt_in, int Mary)
  {
    nt=nt_in;
    set_Gray_QAM(nt,Mary);
  };

  void ND_UQAM::set_Gray_QAM(int nt_in, int Mary)
  {
    nt=nt_in;
    ivec Mary_temp(nt);
    for (int i=0; i<nt; i++) {
      Mary_temp(i)=Mary;
    }
    set_Gray_QAM(nt,Mary_temp);
  };

  void ND_UQAM::set_Gray_QAM(int nt_in, ivec Mary)
  {
    nt=nt_in;
    it_assert(length(Mary)==nt,"ND_UQAM::set_Gray_QAM() Mary has wrong length");
    k.set_size(nt);
    M=to_svec(Mary);
    L.set_size(nt);
    bitmap.set_size(nt);
    symbols.set_size(nt);
    bits2symbols.set_size(nt);    

    for (int i=0; i<nt; i++) {
      k(i) = round_i(log2(double(M(i))));      
      it_assert( ((k(i)>0) && ((1<<k(i))==M(i))),"ND_UQAM::set_Gray_QAM(): M is not a power of 2.");

      L(i) = round_i(std::sqrt((double)M(i)));
      it_assert(L(i)*L(i)== M(i),"ND_UQAM: constellation M must be square");
      
      symbols(i).set_size(M(i)+1);
      bitmap(i).set_size(M(i),k(i));
      bits2symbols(i).set_size(M(i));
      double average_energy = double(M(i)-1)*2.0/3.0;
      double scaling_factor = std::sqrt(average_energy);
      bmat gray_code = graycode(levels2bits(L(i)));
      
      for (int j1=0; j1<L(i); j1++) {
	for (int j2=0; j2<L(i); j2++) {
	  symbols(i)(j1*L(i)+j2) = std::complex<double>( ((L(i)-1)-j2*2.0)/scaling_factor,
							 ((L(i)-1)-j1*2.0)/scaling_factor );
	  bitmap(i).set_row(j1*L(i)+j2, concat(gray_code.get_row(j1), 
					       gray_code.get_row(j2)));
	  bits2symbols(i)( bin2dec(bitmap(i).get_row(j1*L(i)+j2)) ) = j1*L(i)+j2;
	}
      }
      
      symbols(i)(M(i))=0.0;  // must end with a zero; only for a trick exploited in update_norm()
    }    
      
  };

  


  // ----------- MIMO with uniform PSK ---------------

  ND_UPSK::ND_UPSK(int nt_in, int Mary)
  {
    nt=nt_in;
    set_Gray_PSK(nt,Mary);
  };

  void ND_UPSK::set_Gray_PSK(int nt_in, int Mary) 
  { 
    nt=nt_in;
    ivec Mary_temp(nt);
    for (int i=0; i<nt; i++) {
      Mary_temp(i)=Mary;
    }
    set_Gray_PSK(nt,Mary_temp);
  };
  
  void ND_UPSK::set_Gray_PSK(int nt_in, ivec Mary) 
  { 
    nt=nt_in;
    it_assert(length(Mary)==nt,"ND_UPSK::set_Gray_PSK() Mary has wrong length");
    k.set_size(nt);
    M=to_svec(Mary);
    bitmap.set_size(nt);
    symbols.set_size(nt);
    bits2symbols.set_size(nt);    
   
    for (int i=0; i<nt; i++) {
      k(i) = round_i(log2(double(M(i))));      
      it_assert( ((k(i)>0) && ((1<<k(i))==M(i))),"ND_UPSK::set_Gray_PSK(): M is not a power of 2.");
       
      symbols(i).set_size(M(i)+1);
      bits2symbols(i).set_size(M(i));
      bitmap(i) = graycode(k(i));

      double delta = 2.0*pi/M(i);
      double epsilon = delta/10000.0;
      
      for (int j=0; j<M(i); j++) {
	std::complex<double> symb = std::complex<double>(std::polar(1.0,delta*j));

	if (std::abs(std::real(symb)) < epsilon) { 
	  symbols(i)(j) = std::complex<double>(0.0,std::imag(symb)); 
	} else if (std::abs(std::imag(symb)) < epsilon) { 
	  symbols(i)(j) = std::complex<double>(std::real(symb),0.0); }
	else { 
	  symbols(i)(j) = symb; 
	}

	bits2symbols(i)(bin2dec(bitmap(i).get_row(j))) = j;
      }      
      
      symbols(i)(M(i))=0.0;  // must end with a zero; only for a trick exploited in update_norm()
    }


  };
  
  


} // namespace itpp
