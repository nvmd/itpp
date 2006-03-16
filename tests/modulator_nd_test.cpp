/*!
 * \file 
 * \brief Vector ("MIMO") modulator classes test program
 * \author Erik G. Larsson
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

using namespace std;
using namespace itpp;

#if defined(HAVE_LAPACK)

int main()
{
  cout << "========================================================" << endl;
  cout << "               Test of ND (MIMO) Modulators             " << endl;
  cout << "========================================================" << endl;
   
  RNG_reset(12345);

  {
    ND_UPAM chan;
    cout << chan.get_llrcalc() << endl;
    double sigma2=0.005;
    int nt=5;
    for (int np=1; np<=3; np++) {
      cout << "================== ND-U" << (1<<np) << "PAM =========== " 
	   << endl;

      chan.set_Gray_PAM(nt,1<<np);   
      bvec b=randb(nt*np);
      QLLRvec LLR_ap = zeros_i(nt*np);
      //       LLR_ap(3)=chan.get_llrcalc().to_qllr(100);
      //       LLR_ap(4)=chan.get_llrcalc().to_qllr(-100);
      cout << b << endl;
      vec x=chan.modulate_bits(b);
      
      QLLRvec LLR(nt*np);
      mat H = randn(nt,nt);
      vec y= H*x;
      y+=sqrt(sigma2)*randn(nt);
      //      cout << y << endl;
      chan.map_demod(LLR_ap, LLR, sigma2, H, y);
      cout << "full channel:" << chan.get_llrcalc().to_double(LLR) << endl;

      vec h =diag(H);
      chan.map_demod(LLR_ap, LLR, sigma2, h, y);
      cout << "diagonal channel: " << chan.get_llrcalc().to_double(LLR) << endl;

      //       int bit=0;
      //       for (int i=0; i<nt; i++) {
      // 	ND_UPAM scalchan(1,1<<np);
      // 	mat H0(1,1);
      // 	H0(0,0)=H(i,i);
      // 	vec Y0(1);
      // 	Y0(0)=y(i);
      // 	QLLRvec llrapr = LLR_ap.mid(bit,chan.get_k()(i));
      // 	QLLRvec llrapost(chan.get_k()(i));
      // 	scalchan.map_demod(llrapr,llrapost, sigma2, H0, Y0);
      // 	bit+=chan.get_k()(i);
      // 	cout << chan.get_llrcalc().to_double(llrapost) << endl;
      //       }

      ivec zhat;
      chan.sphere_decoding(y,H, 0.01, 10000, 2.0, zhat);
      cout << zhat << endl;
    }
  }

  {
    ND_UQAM chan;
    double sigma2=0.005;
    int nt=3;
    for (int np=1; np<=3; np++) {
      cout << "================== ND-U" << ((1<<(2*np))) << "QAM =========== " 
	   << endl;

      chan.set_Gray_QAM(nt,(1<<(2*np)));   
      bvec b=randb(nt*np*2);
      cout << b << endl;
      cvec x=chan.modulate_bits(b);

      QLLRvec LLR_ap = zeros_i(2*nt*np);
      QLLRvec LLR(2*nt*np);
      cmat H = randn_c(nt,nt);
      cvec y= H*x;
      y+=sqrt(sigma2)*randn_c(nt);
      chan.map_demod(LLR_ap, LLR, sigma2, H, y);
      cout << "full channel:" << chan.get_llrcalc().to_double(LLR) << endl;

      cvec h =diag(H);
      chan.map_demod(LLR_ap, LLR, sigma2, h, y);
      cout << "diagonal channel: " << chan.get_llrcalc().to_double(LLR) << endl;
    }
  }

  {
    ND_UPSK chan;
    double sigma2=0.005;
    int nt=3;
    for (int np=1; np<=3; np++) {
      cout << "================== ND-U" << ((1<<(2*np))) << "PSK =========== " 
	   << endl;

      chan.set_Gray_PSK(nt,(1<<(2*np)));   
      bvec b=randb(nt*np*2);
      cout << b << endl;
      cvec x=chan.modulate_bits(b);

      QLLRvec LLR_ap = zeros_i(2*nt*np);
      QLLRvec LLR(2*nt*np);
      cmat H = randn_c(nt,nt);
      cvec y= H*x;
      y+=sqrt(sigma2)*randn_c(nt);
      chan.map_demod(LLR_ap, LLR, sigma2, H, y);
      cout << "full channel:" << chan.get_llrcalc().to_double(LLR) << endl;

      cvec h =diag(H);
      chan.map_demod(LLR_ap, LLR, sigma2, h, y);
      cout << "diagonal channel: " << chan.get_llrcalc().to_double(LLR) << endl;
    }
  }
}

#else

int main() { 
  cerr << "Error: LAPACK library is needed to run this test program" << endl;
  return 1;
}

#endif
