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

#ifndef ITFILE_TEST_FILE

int main() {
  cout << "ITFILE_TEST_FILE not defined. Test skipped." << endl;
}

#else

int main()
{
  string c = "abcdefgh";
  string d;
  int x = 1234567890;
  int y;

  it_file_base::data_header dh;

  it_file ff;
  ff.open(string(ITFILE_TEST_FILE));
#ifdef SAVE_DATA  
  ff << Name("x") << x;
  ff << Name("c") << c;
#endif
  ff >> Name("x") >> y;
  ff >> Name("c") >> d;
  ff.close();

  ff.open(string(ITFILE_TEST_FILE));
  ff.read_data_header(dh);
  ff.close();

  if (dh.endianity == 0)
    cout << "Little endian data:" << endl;
  else if (dh.endianity == 1)
    cout << "Big endian data:" << endl;
  else
    it_error("Endianness error");
  
  cout << "x = " << x << endl;
  cout << "y = " << y << endl;
  cout << "c = " << c << endl;
  cout << "d = " << d << endl;

  return 0;
}

#endif
