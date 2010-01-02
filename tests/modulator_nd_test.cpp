/*!
 * \file
 * \brief Vector ("MIMO") modulator classes test program
 * \author Erik G. Larsson and Adam Piatyszek
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
#include <iomanip>

using namespace std;
using namespace itpp;


int main()
{
  cout << "========================================================" << endl;
  cout << "               Test of ND (MIMO) Modulators             " << endl;
  cout << "========================================================" << endl;

  cout.setf(ios::fixed);
  cout.precision(2);

  RNG_reset(12345);
  double sigma2 = 0.005;
  double sigma = std::sqrt(sigma2);

  {
    ND_UPAM chan;
    int nt = 5;
    for (int np = 1; np <= 3; np++) {
      cout << "================== ND-U" << (1 << np) << "PAM ==================\n";

      chan.set_M(nt, 1 << np);
      cout << chan << endl;
      bvec b = randb(nt * np);
      cout << b << endl;
      vec x = chan.modulate_bits(b);
      mat H = randn(nt, nt);
      vec y = H * x + sigma * randn(nt);

      QLLRvec LLR_ap = zeros_i(nt * np);
      QLLRvec LLR;

      chan.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR);
      cout << "full channel     : " << chan.get_llrcalc().to_double(LLR) << endl;
      chan.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR, ND_UPAM::FULL_ENUM_LOGMAP);
      cout << "                   " << chan.get_llrcalc().to_double(LLR) << endl;

      chan.demodulate_soft_bits(y, diag(H), sigma2, LLR_ap, LLR);
      cout << "diagonal channel : " << chan.get_llrcalc().to_double(LLR) << endl;

      chan.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR, ND_UPAM::ZF_LOGMAP);
      cout << "zero-forcing     : " << chan.get_llrcalc().to_double(LLR) << endl;

      ivec zhat;
      chan.sphere_decoding(y, H, 0.01, 10000, 2.0, zhat);
      cout << zhat << endl;
    }
  }

  {
    ND_UQAM chan;
    int nt = 3;
    for (int np = 1; np <= 3; np++) {
      cout << "================== ND-U" << ((1 << (2*np))) << "QAM ==================\n";

      chan.set_M(nt, (1 << (2*np)));
      cout << chan << endl;
      bvec b = randb(nt * np * 2);
      cout << b << endl;
      cvec x = chan.modulate_bits(b);
      cmat H = randn_c(nt, nt);
      cvec y = H * x + sigma * randn_c(nt);
      QLLRvec LLR_ap = zeros_i(2 * nt * np);
      QLLRvec LLR(2*nt*np);

      chan.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR);
      cout << "full channel     : " << chan.get_llrcalc().to_double(LLR) << endl;
      chan.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR, ND_UQAM::FULL_ENUM_LOGMAP);
      cout << "                 : " << chan.get_llrcalc().to_double(LLR) << endl;

      chan.demodulate_soft_bits(y, diag(H), sigma2, LLR_ap, LLR);
      cout << "diagonal channel : " << chan.get_llrcalc().to_double(LLR) << endl;

      chan.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR, ND_UPAM::ZF_LOGMAP);
      cout << "zero-forcing     : " << chan.get_llrcalc().to_double(LLR) << endl;
    }
  }

  {
    ND_UPSK chan;
    int nt = 3;
    for (int np = 1; np <= 3; np++) {
      cout << "================== ND-U" << ((1 << (2*np))) << "PSK ==================\n";

      chan.set_M(nt, (1 << (2*np)));
      cout << chan << endl;
      bvec b = randb(nt * np * 2);
      cout << b << endl;
      cvec x = chan.modulate_bits(b);
      cmat H = randn_c(nt, nt);
      cvec y = H * x + sigma * randn_c(nt);
      QLLRvec LLR_ap = zeros_i(2 * nt * np);
      QLLRvec LLR(2*nt*np);

      chan.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR);
      cout << "full channel     : " << chan.get_llrcalc().to_double(LLR) << endl;
      chan.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR, ND_UQAM::FULL_ENUM_LOGMAP);
      cout << "                 : " << chan.get_llrcalc().to_double(LLR) << endl;

      chan.demodulate_soft_bits(y, diag(H), sigma2, LLR_ap, LLR);
      cout << "diagonal channel : " << chan.get_llrcalc().to_double(LLR) << endl;

      chan.demodulate_soft_bits(y, H, sigma2, LLR_ap, LLR, ND_UPAM::ZF_LOGMAP);
      cout << "zero-forcing     : " << chan.get_llrcalc().to_double(LLR) << endl;
    }
  }

  return 0;
}
