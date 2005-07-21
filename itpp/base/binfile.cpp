/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2003 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Binary file formats implementations
  \author Tony Ottosson and Thomas Eriksson

  $Revision$

  $Date$
*/

#include <itpp/base/binfile.h>
#include <itpp/base/machdep.h>
#include <fstream>
#include <sys/stat.h>
#include <string>

using std::string;
using std::ofstream;
using std::ifstream;
using std::fstream;
using std::ios;

namespace itpp { 

  bool exist(const std::string &name)
  {
    struct stat st;

    return stat(name.c_str(), &st) == 0;
  }

  bfstream_base::bfstream_base(endian e)
  {
    endianity = e;
    native_endianity = (big_endian((short)0x1234) == (short)0x1234)
      ? b_endian : l_endian;
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
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity) {
      put(c[0]);
      put(c[1]);
    } else {
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bofstream& bofstream::operator <<(unsigned short a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 2);
    else {
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bofstream& bofstream::operator <<(float a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 4);
    else {
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bofstream& bofstream::operator <<(double a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 8);
    else {
      put(c[7]);
      put(c[6]);
      put(c[5]);
      put(c[4]);
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

//   bofstream& bofstream::operator <<(long double a)
//   {
//     char *c=reinterpret_cast<char *>(&a);

//     if (endianity == native_endianity)
//       write(c, 10);
//     else {
//       put(c[9]);
//       put(c[8]);
//       put(c[7]);
//       put(c[6]);
//       put(c[5]);
//       put(c[4]);
//       put(c[3]);
//       put(c[2]);
//       put(c[1]);
//       put(c[0]);
//     }
//     return *this;
//   }

  bofstream& bofstream::operator <<(int a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 4);
    else {
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bofstream& bofstream::operator <<(unsigned a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 4);
    else {
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bofstream& bofstream::operator <<(long a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 4);
    else {
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bofstream& bofstream::operator <<(unsigned long a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 4);
    else {
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
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
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 4);
    else {
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bifstream& bifstream::operator >>(unsigned int &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 4);
    else {
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bifstream& bifstream::operator >>(short int &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 2);
    else {
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bifstream& bifstream::operator >>(unsigned short int &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 2);
    else {
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bifstream& bifstream::operator >>(float &a)
  {
    char *c=reinterpret_cast<char *>(&a);

	if (endianity == native_endianity) {
      read(c, 4);
	} else {
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bifstream& bifstream::operator >>(double &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 8);
    else {
      get(c[7]);
      get(c[6]);
      get(c[5]);
      get(c[4]);
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

//   bifstream& bifstream::operator >>(long double &a)
//   {
//     char *c=reinterpret_cast<char *>(&a);

//     if (endianity == native_endianity)
//       read(c, 10);
//     else {
//       get(c[9]);
//       get(c[8]);
//       get(c[7]);
//       get(c[6]);
//       get(c[5]);
//       get(c[4]);
//       get(c[3]);
//       get(c[2]);
//       get(c[1]);
//       get(c[0]);
//     }
//     return *this;
//   }

  bifstream& bifstream::operator >>(long int &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 4);
    else {
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bifstream& bifstream::operator >>(unsigned long int &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 4);
    else {
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
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

  void bfstream::open(const std::string &name, bool trnc, endian e) //CC fix trunc -> trnc
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
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 4);
    else {
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator <<(unsigned int a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 4);
    else {
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator <<(short a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity) {
      put(c[0]);
      put(c[1]);
    } else {
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator <<(unsigned short a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 2);
    else {
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator <<(float a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 4);
    else {
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator <<(double a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 8);
    else {
      put(c[7]);
      put(c[6]);
      put(c[5]);
      put(c[4]);
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

//   bfstream& bfstream::operator <<(long double a)
//   {
//     char *c=reinterpret_cast<char *>(&a);

//     if (endianity == native_endianity)
//       write(c, 10);
//     else {
//       put(c[9]);
//       put(c[8]);
//       put(c[7]);
//       put(c[6]);
//       put(c[5]);
//       put(c[4]);
//       put(c[3]);
//       put(c[2]);
//       put(c[1]);
//       put(c[0]);
//     }
//     return *this;
//   }

  bfstream& bfstream::operator <<(long int a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 4);
    else {
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator <<(unsigned long int a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      write(c, 4);
    else {
      put(c[3]);
      put(c[2]);
      put(c[1]);
      put(c[0]);
    }
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
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 4);
    else {
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator >>(unsigned int &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 4);
    else {
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator >>(short int &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 2);
    else {
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator >>(unsigned short int &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 2);
    else {
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator >>(float &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 4);
    else {
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator >>(double &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 8);
    else {
      get(c[7]);
      get(c[6]);
      get(c[5]);
      get(c[4]);
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

//   bfstream& bfstream::operator >>(long double &a)
//   {
//     char *c=reinterpret_cast<char *>(&a);

//     if (endianity == native_endianity)
//       read(c, 10);
//     else {
//       get(c[9]);
//       get(c[8]);
//       get(c[7]);
//       get(c[6]);
//       get(c[5]);
//       get(c[4]);
//       get(c[3]);
//       get(c[2]);
//       get(c[1]);
//       get(c[0]);
//     }
//     return *this;
//   }

  bfstream& bfstream::operator >>(long int &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 4);
    else {
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
    return *this;
  }

  bfstream& bfstream::operator >>(unsigned long int &a)
  {
    char *c=reinterpret_cast<char *>(&a);

    if (endianity == native_endianity)
      read(c, 4);
    else {
      get(c[3]);
      get(c[2]);
      get(c[1]);
      get(c[0]);
    }
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

} //namespace itpp
