/*!
 * \file 
 * \brief Binary file formats implementations
 * \author Tony Ottosson, Thomas Eriksson and Adam Piatyszek
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

#include <itpp/base/binfile.h>
#include <itpp/base/it_endian.h>

using std::ofstream;
using std::ifstream;
using std::fstream;
using std::ios;


namespace itpp { 

  bool exist(const std::string &name)
  {
    bool file_exists = false;
    ifstream file(name.c_str(), ios::in);
    if (file.is_open()) {
      file_exists = true;
    }
    file.close();
    return file_exists;
  }

  bfstream_base::bfstream_base(endian e)
  {
    endianity = e;
    native_endianity = (check_big_endianness() ? b_endian : l_endian);
  }

  //-----------------------------------------------------------------------
  //	bofstream
  //-----------------------------------------------------------------------

  bofstream::bofstream(const std::string &name, endian e)
    : bfstream_base(e), ofstream(name.c_str(), ios::out | ios::binary)
  {
  }

  bofstream::bofstream()
    : bfstream_base(), ofstream()
  {
  }

  void bofstream::open(const std::string &name, endian e)
  {
    endianity = e;
    ofstream::open(name.c_str(), ios::out | ios::binary);
  }

  bofstream& bofstream::operator <<(char a)
  {
    put(a);
    return *this;
  }

  bofstream& bofstream::operator <<(const class bin &a)
  {
    put(a.value());
    return *this;
  }

  bofstream& bofstream::operator <<(short a)
  {
    write_endian<short, bofstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bofstream& bofstream::operator <<(unsigned short a)
  {
    write_endian<unsigned short, bofstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bofstream& bofstream::operator <<(float a)
  {
    write_endian<float, bofstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bofstream& bofstream::operator <<(double a)
  {
    write_endian<double, bofstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bofstream& bofstream::operator <<(int a)
  {
    write_endian<int, bofstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bofstream& bofstream::operator <<(unsigned int a)
  {
    write_endian<unsigned int, bofstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bofstream& bofstream::operator <<(long int a)
  {
    write_endian<long int, bofstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bofstream& bofstream::operator <<(unsigned long int a)
  {
    write_endian<unsigned long int, bofstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bofstream& bofstream::operator<<(const char *a)
  {
    write(a, strlen(a)+1);
    return *this;
  }

  bofstream& bofstream::operator<<(const std::string &a)
  {
    write(a.c_str(), a.size()+1);
    return *this;
  }

  //-----------------------------------------------------------------------
  //	bifstream
  //-----------------------------------------------------------------------

  bifstream::bifstream(const std::string &name, endian e)
    : bfstream_base(e), ifstream(name.c_str(), ios::in | ios::binary )
  {
  }

  bifstream::bifstream()
    : bfstream_base(), ifstream()
  {
  }

  void bifstream::open(const std::string &name, endian e)
  {
    endianity = e;
    ifstream::open(name.c_str(),ios::in | ios::binary );
  }

  long bifstream::length()  //in bytes
  {
    std::streampos pos1,len;
    pos1=tellg();
    seekg(0,ios::end);
    len=tellg();
    seekg(pos1);
    return len;
  }

  bifstream& bifstream::operator >>(char &a)
  {
    get(a);
    return *this;
  }

  bifstream& bifstream::operator >>(class bin &a)
  {
    char temp;
    get(temp);
    a=temp;
    return *this;
  }

  bifstream& bifstream::operator >>(int &a)
  {
    a = read_endian<int, bifstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bifstream& bifstream::operator >>(unsigned int &a)
  {
    a = read_endian<unsigned int, bifstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bifstream& bifstream::operator >>(short int &a)
  {
    a = read_endian<short int, bifstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bifstream& bifstream::operator >>(unsigned short int &a)
  {
    a = read_endian<unsigned short int, bifstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bifstream& bifstream::operator >>(float &a)
  {
    a = read_endian<float, bifstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bifstream& bifstream::operator >>(double &a)
  {
    a = read_endian<double, bifstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bifstream& bifstream::operator >>(long int &a)
  {
    a = read_endian<long int, bifstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bifstream& bifstream::operator >>(unsigned long int &a)
  {
    a = read_endian<unsigned long int, bifstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bifstream& bifstream::operator>>(char *a)
  {
    getline(a, '\0');
    return *this;
  }

  bifstream& bifstream::operator>>(std::string &a)
  {
    std::getline(*this, a, '\0');
    return *this;
  }

  //-----------------------------------------------------------------------
  //	bfstream
  //-----------------------------------------------------------------------

  bfstream::bfstream(const std::string &name, endian e)
    : bfstream_base(e), fstream(name.c_str(), ios::in | ios::out | ios::binary)
  {
  }

  bfstream::bfstream()
    : bfstream_base(), fstream()
  {
  }

  void bfstream::open(const std::string &name, bool trnc, endian e)
  {
    endianity = e;

    if (trnc)
      fstream::open(name.c_str(), ios::in | ios::out | ios::binary | ios::trunc);
    else
      fstream::open(name.c_str(), ios::in | ios::out | ios::binary);
  }

  void bfstream::open_readonly(const std::string &name, endian e)
  {
    endianity = e;
    fstream::open(name.c_str(), ios::in | ios::binary);
  }

  long bfstream::length()  //in bytes
  {
    std::streampos pos1,len;
    pos1=tellg();
    seekg(0,ios::end);
    len=tellg();
    seekg(pos1);
    return len;
  }

  bfstream& bfstream::operator <<(char a)
  {
    put(a);
    return *this;
  }

  bfstream& bfstream::operator <<(const class bin &a)
  {
    put(a.value());
    return *this;
  }

  bfstream& bfstream::operator <<(int a)
  {
    write_endian<int, bfstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator <<(unsigned int a)
  {
    write_endian<unsigned int, bfstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator <<(short a)
  {
    write_endian<short, bfstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator <<(unsigned short a)
  {
    write_endian<unsigned short, bfstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator <<(float a)
  {
    write_endian<float, bfstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator <<(double a)
  {
    write_endian<double, bfstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator <<(long int a)
  {
    write_endian<long int, bfstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator <<(unsigned long int a)
  {
    write_endian<unsigned long int, bfstream&>(*this, a, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator<<(const char *a)
  {
    write(a, strlen(a)+1);
    return *this;
  }

  bfstream& bfstream::operator<<(const std::string &a)
  {
    write(a.c_str(), a.size()+1);
    return *this;
  }

  bfstream& bfstream::operator >>(char &a)
  {
    get(a);
    return *this;
  }

  bfstream& bfstream::operator >>(class bin &a)
  {
    char temp;
    get(temp);
    a=temp;
    return *this;
  }

  bfstream& bfstream::operator >>(int &a)
  {
    a = read_endian<int, bfstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator >>(unsigned int &a)
  {
    a = read_endian<unsigned int, bfstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator >>(short int &a)
  {
    a = read_endian<short int, bfstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator >>(unsigned short int &a)
  {
    a = read_endian<unsigned short int, bfstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator >>(float &a)
  {
    a = read_endian<float, bfstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator >>(double &a)
  {
    a = read_endian<double, bfstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator >>(long int &a)
  {
    a = read_endian<long int, bfstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator >>(unsigned long int &a)
  {
    a = read_endian<unsigned long int, bfstream&>(*this, endianity != native_endianity);
    return *this;
  }

  bfstream& bfstream::operator>>(char *a)
  {
    getline(a, '\0');
    return *this;
  }

  bfstream& bfstream::operator>>(std::string &a)
  {
    std::getline(*this, a, '\0');
    return *this;
  }

} // namespace itpp
