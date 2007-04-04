/*!
* \file 
* \brief IT file endianness test program
* \author Adam Piatyszek
*
* $Date$
* $Revision$
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

#include <itpp/itbase.h>

using namespace itpp;
using namespace std;

// To run extensive tests uncomment the following definition
// #define EXTENSIVE_TESTS

#ifndef ITFILE_TEST_FILE

int main() {
  cout << "ITFILE_TEST_FILE not defined. Test skipped." << endl;
}

#else

int main()
{
  char c, c_ref = 'c';
  bin b, b_ref = 0;
  short s, s_ref = 1234;
  int i, i_ref = -1234567890;
  float f, f_ref = -12345.6;
  double d, d_ref = 2.1e-8;
  complex<double> cd, cd_ref = std::complex<double>(1.0, -1.5);
  string st, st_ref = "abcdefghij 0123456789";

  bvec bv, bv_ref = "0 1 0 1 1";
  ivec iv, iv_ref = "2 5 -400 2 -10";
  vec v, v_ref = "1e-9 0.2 0.7 -1.0 0.0";
  cvec cv, cv_ref = "(0,2) (1.5,7.2)";

  bmat bm, bm_ref = "0 1 0; 1 1 1";
  imat im, im_ref = "2 5; 4 -10; 0 3";
  mat m, m_ref = "1e-9 0.2 0.7; -1.0 0.0 3e10";
  cmat cm, cm_ref = "(0,2) (-1.5,7.2); (1.1,2) (7,-4e-5)";

  Array<bvec> abv, abv_ref = "{[0 1] [0 0 0] [1 0] [1 1 0]}";
  Array<ivec> aiv, aiv_ref = "{[1 2 3] [4 5 6] [7 8 9]}";
  Array<vec> av, av_ref = "{[1 2e4 -3.0] [-0.5 6e-9] [7e-3]}";
  Array<cvec> acv, acv_ref = "{[(0,2) (0.5,-0.5)] [(2,1)] [(0,0) (-4.6,2)]}";

  Array<bmat> abm, abm_ref = "{[0 1 0; 0 0 1] [1 1; 0 0; 1 0] [1; 1; 0]}";
  Array<imat> aim, aim_ref = "{[0 2 3; 0 -1 9; 2 3 -1] [1 10 100] [0 4; -2 3]}";
  Array<mat> am, am_ref = "{[0.5 2e7; 0.5 -0.5] [1e-4 3 4; 0.1 0.2 .3]}";
  Array<cmat> acm, acm_ref = "{[(0,2) (0.5,-0.5); (2,1) (0,0)] "
    "[(1.1,2) (7,-4e-5); (0,2) (1.5,7.2)]}";

  it_file ff;
  ff.open(string(ITFILE_TEST_FILE));
#ifdef SAVE_DATA  
  ff << Name("c") << c_ref;
  ff << Name("b") << b_ref;
  ff << Name("s") << s_ref;
  ff << Name("i") << i_ref;
  ff << Name("f") << f_ref;
  ff << Name("d") << d_ref;
  ff << Name("cd") << cd_ref;
  ff << Name("st") << st_ref;
  ff << Name("bv") << bv_ref;
  ff << Name("iv") << iv_ref;
  ff << Name("v") << v_ref;
  ff << Name("cv") << cv_ref;
  ff << Name("bm") << bm_ref;
  ff << Name("im") << im_ref;
  ff << Name("m") << m_ref;
  ff << Name("cm") << cm_ref;
  ff << Name("abv") << abv_ref;
  ff << Name("aiv") << aiv_ref;
  ff << Name("av") << av_ref;
  ff << Name("acv") << acv_ref;
  ff << Name("abm") << abm_ref;
  ff << Name("aim") << aim_ref;
  ff << Name("am") << am_ref;
  ff << Name("acm") << acm_ref;
#endif
  ff >> Name("c") >> c;
  ff >> Name("b") >> b;
  ff >> Name("s") >> s;
  ff >> Name("i") >> i;
  ff >> Name("f") >> f;
  ff >> Name("d") >> d;
  ff >> Name("cd") >> cd;
  ff >> Name("st") >> st;
  ff >> Name("bv") >> bv;
  ff >> Name("iv") >> iv;
  ff >> Name("v") >> v;
  ff >> Name("cv") >> cv;
  ff >> Name("bm") >> bm;
  ff >> Name("im") >> im;
  ff >> Name("m") >> m;
  ff >> Name("cm") >> cm;
  ff >> Name("abv") >> abv;
  ff >> Name("aiv") >> aiv;
  ff >> Name("av") >> av;
  ff >> Name("acv") >> acv;
  ff >> Name("abm") >> abm;
  ff >> Name("aim") >> aim;
  ff >> Name("am") >> am;
  ff >> Name("acm") >> acm;
  ff.close();

  cout << "char    : '" << c << "'"<< endl
       << "          '" << c_ref << "'" << endl
       << "bin     : " << b << endl
       << "          " << b_ref << endl
       << "short   : " << s << endl
       << "          " << s_ref << endl
       << "int     : " << i << endl
       << "          " << i_ref << endl
       << "float   : " << f << endl
       << "          " << f_ref << endl
       << "double  : " << d << endl
       << "          " << d_ref << endl
       << "complex : " << cd << endl
       << "          " << cd_ref << endl
       << "string  : \"" << st << "\"" << endl
       << "          \"" << st_ref << "\"" << endl << endl;

  cout << "bvec    : " << bv << endl
       << "          " << bv_ref << endl
       << "ivec    : " << iv << endl
       << "          " << iv_ref << endl
       << "vec     : " << v << endl
       << "          " << v_ref << endl
       << "cvec    : " << cv << endl
       << "          " << cv_ref << endl << endl;

  cout << "bmat    :\n" << bm << endl << bm_ref << endl
       << "imat    :\n" << im << endl << im_ref << endl
       << "mat     :\n" << m << endl << m_ref << endl
       << "cmat    :\n" << cm << endl << cm_ref << endl << endl;

  cout << "Array<bvec> :\n" << abv << endl << abv_ref << endl
       << "Array<ivec> :\n" << aiv << endl << aiv_ref << endl
       << "Array<vec>  :\n" << av << endl << av_ref << endl
       << "Array<cvec> :\n" << acv << endl << acv_ref << endl << endl;

  cout << "Array<bmat> :\n" << abm << endl << abm_ref << endl
       << "Array<imat> :\n" << aim << endl << aim_ref << endl
       << "Array<mat> :\n" << am << endl << am_ref << endl
       << "Array<cmat> :\n" << acm << endl << acm_ref << endl << endl;

#ifdef EXTENSIVE_TESTS
  ivec iv0 = "0 0";
  ivec iv1 = ones_i(100);
  ivec iv2 = "2 2 2 2";
  ivec iv3 = "3";

  it_file ff1("itfile_test_extensive.it", true);
  ff1 << Name("iv0") << iv0 << flush;
  ff1 << Name("iv1") << iv1 << flush;
  ff1 << Name("iv2") << iv2 << flush;
  ff1.remove("iv1");
  ff1 << Name("iv1") << ivec("1");
  ff1 << Name("iv3") << iv3 << flush;
  ff1 << Name("iv4") << iv3 << flush;
  ff1.close();

  std::string name, type;
  int size;
  it_ifile ff2("itfile_test_extensive.it");
  int n = 0;
  while (ff2.seek(n++)) {
    ff2.info(name, type, size);
    cout << "Name = " << name << "   Type = " << type << "   Size = " << size 
	 << endl;
    ff2 >> iv1;
    cout << iv1 << endl;
  }
  ff2.close();
#endif

  return 0;
}

#endif
