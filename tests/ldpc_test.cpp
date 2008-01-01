/*!
 * \file
 * \brief LDPC class test program
 * \author Erik G. Larsson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2008  (see AUTHORS file for a list of contributors)
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

using std::cout;
using std::endl;
using namespace itpp;

int main()
{
  LDPC_Parity_Regular H;
  H.generate(200, 3, 6, "rand", "100 6");
  int girth = H.cycle_removal_MGW(6);
  cout << "girth=" << girth << endl;
  LDPC_Generator_Systematic G;
  G.construct(&H);
  LDPC_Code C(&H, &G);
  C.save_code("ldpc_test.codec");
  LDPC_Code C1("ldpc_test.codec", &G);
  cout << C << endl;
  bvec bitsin=randb(C.get_nvar()-C.get_ncheck());
  bvec bitsout;
  C.encode(bitsin,bitsout);
  it_assert(C.syndrome_check(bitsout),"syndrome check failed");

  double EbN0db=1.5;
  double N0 = pow(10.0,-EbN0db/10.0) / C.get_rate();
  double sigma=sqrt(N0/2.0);
  vec x = 1.0+sigma*randn(C.get_nvar());
  QLLRvec LLRin = C.get_llrcalc().to_qllr(2.0*x/(N0/2.0));
  QLLRvec LLRout(C.get_nvar());
  C.bp_decode(LLRin,LLRout);
  cout << LLRout.left(25) << endl;
}
