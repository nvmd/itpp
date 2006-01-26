/*!
 * \file 
 * \brief Implementation of audio Audio classes and functions
 * \author Tobias Ringstrom
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

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#include <itpp/srccode/audiofile.h>
#include <itpp/base/machdep.h>
#include <iostream>


using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::ios;

namespace itpp {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#define SND_MAGIC    0x2e736e64

  inline static short double_to_short(double x)
  {
    if (x >= 32767.0)
      return 32767;
    else if (x <= -32768.0)
      return -32768;
    else
      return round_i(x);
  }

  inline static signed char double_to_char(double x)
  {
    if (x >= 127.0)
      return 127;
    else if (x <= -128.0)
      return -128;
    else
      return round_i(x);
  }

#define MK_RW_FUNS(type,namestub)					\
  inline static type read_##namestub(istream &s)			\
  { type v; s.read(reinterpret_cast<char *>(&v), sizeof(type)); return v; } \
    inline static void write_##namestub(ostream &s, type v)		\
    { s.write(reinterpret_cast<char *>(&v), sizeof(type)); }

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

  MK_RW_FUNS(signed char,char)
    MK_RW_FUNS(short,short)
    MK_RW_FUNS(int,int)
    MK_RW_FUNS(unsigned,unsigned)
    MK_RW_FUNS(float,float)
    MK_RW_FUNS(double,double)

    bool raw16le_read(const char *fname, vec &v)
  {
    ifstream file(fname, ios::in | ios::binary);
    int n, i;

#ifdef HAVE_SYS_STAT_H
    struct stat st;
    if (stat(fname, &st) == -1)
      return false;
#endif

    n = st.st_size / 2; // short vs. byte
    v.set_size(n, false);
    for (i=0; i<n; i++)
      v(i) = little_endian(read_short(file)) / 32768.0;
    if (!file)
      return false;
	
    return true;
  }

  bool raw16le_read(const char *fname, vec &v, int beg, int len)
  {
    ifstream file(fname, ios::in | ios::binary);
    int i;

    it_assert1(len >= 0, "raw16_read()");
    v.set_size(len, false);
    file.seekg(2 * beg);
    for (i=0; i<len; i++)
      v(i) = little_endian(read_short(file)) / 32768.0;
    if (!file)
      return false;
	
    return true;
  }

  bool raw16le_write(const char *fname, const vec &v, bool append)
  {
    ofstream file(fname, (append ? ios::app | ios::ate : ios::out | ios::trunc) | ios::binary);
    int i;

    for (i=0; i<v.size(); i++)
      write_short(file, little_endian(double_to_short(v(i) * 32768.0)));
    if (!file)
      return false;
	
    return true;
  }

  bool raw16be_read(const char *fname, vec &v)
  {
    ifstream file(fname, ios::in | ios::binary);
    struct stat st;
    int n, i;

    if (stat(fname, &st) == -1)
      return false;

    n = st.st_size / 2; // short vs. byte
    v.set_size(n, false);
    for (i=0; i<n; i++)
      v(i) = big_endian(read_short(file)) / 32768.0;
    if (!file)
      return false;
	
    return true;
  }

  bool raw16be_read(const char *fname, vec &v, int beg, int len)
  {
    ifstream file(fname, ios::in | ios::binary);
    int i;

    it_assert1(len >= 0, "raw16_read()");
    v.set_size(len, false);
    file.seekg(2 * beg);
    for (i=0; i<len; i++)
      v(i) = big_endian(read_short(file)) / 32768.0;
    if (!file)
      return false;
	
    return true;
  }

  bool raw16be_write(const char *fname, const vec &v, bool append)
  {
    ofstream file(fname, (append ? ios::app | ios::ate : ios::out | ios::trunc) | ios::binary);
    int i;

    for (i=0; i<v.size(); i++)
      write_short(file, big_endian(double_to_short(v(i) * 32768.0)));
    if (!file)
      return false;
	
    return true;
  }

  //////////////////////////////////////////////////
  //
  // Audio_File
  //
  //////////////////////////////////////////////////
  Audio_File::Audio_File()
  {
    is_valid = false;
  }

  //////////////////////////////////////////////////
  //
  // SND_Format
  //
  //////////////////////////////////////////////////
  int SND_Format::sample_size() const
  {
    switch (header.encoding) {
    case enc_mulaw8 :   return 1;
    case enc_alaw8 :    return 1;
    case enc_linear8 :  return 1;
    case enc_linear16 : return 2;
    case enc_linear24 : return 3;
    case enc_linear32 : return 4;
    case enc_float :    return 4;
    case enc_double :   return 8;
    }
    return 0;
  }

  bool SND_Format::read_header(istream &f)
  {
    f.seekg(0);
    header.magic = big_endian(read_unsigned(f));
    header.hdr_size = big_endian(read_unsigned(f));
    header.data_size = big_endian(read_unsigned(f));
    header.encoding = big_endian(read_unsigned(f));
    header.sample_rate = big_endian(read_unsigned(f));
    header.channels = big_endian(read_unsigned(f));
    f.read(header.info, SND_INFO_LEN);
    if (!f || header.magic != SND_MAGIC) {
      std::cerr << header.magic << " != " << SND_MAGIC << std::endl;
      it_warning("SND_Format::read_header(): This is not a .snd file!");
      return false;
    }
    f.seekg(header.hdr_size);

    return f.good();
  }

  bool SND_Format::write_header(ostream &f)
  {
    f.seekp(0);
    header.magic = SND_MAGIC;
    header.hdr_size = sizeof(header);
    memset(header.info, 0, SND_INFO_LEN);

    write_unsigned(f, big_endian(header.magic));
    write_unsigned(f, big_endian(header.hdr_size));
    write_unsigned(f, big_endian(header.data_size));
    write_unsigned(f, big_endian(header.encoding));
    write_unsigned(f, big_endian(header.sample_rate));
    write_unsigned(f, big_endian(header.channels));
    f.write(reinterpret_cast<char *>(&header.info), SND_INFO_LEN);

    return f.good();
  }

  //////////////////////////////////////////////////
  //
  // SND_In_File
  //
  //////////////////////////////////////////////////

  SND_In_File::SND_In_File()
  {
  }

  SND_In_File::SND_In_File(const char *fname)
  {
    open(fname);
  }

  bool SND_In_File::open(const char *fname)
  {
    if (file.is_open())
      close();
    file.clear();
    is_valid = false;
    file.open(fname, ios::in | ios::binary);
    if (!file)
      return false;
    if (!read_header(file)) {
      file.close();
      return false;
    }

    is_valid = true;
    return true;
  }

  void SND_In_File::close()
  {
    file.close();
    is_valid = false;
  }
    
  bool SND_In_File::seek_read(int pos)
  {
    if (pos < 0)
      file.seekg(0, ios::end);
    else
      file.seekg(header.hdr_size + header.channels * sample_size() * pos);
    return true;
  }

  int SND_In_File::tell_read()
  {
    if (!good())
      return -1;
    
    return ( static_cast<int>(file.tellg()) - sizeof(header) ) / ( header.channels * sample_size() );
  }

  bool SND_In_File::read(vec &v)
  {
    if (!good())
      return false;
    
    int i, n;
    
    n = samples();
    v.set_size(n, false);
    seek_read(0);
    
    switch (header.encoding) {
    case enc_linear8 :
      for (i=0; i<n; i++)
	v(i) = read_char(file) / 128.0;
      break;
    case enc_linear16 :
      for (i=0; i<n; i++)
	v(i) = big_endian(read_short(file)) / 32768.0;
      break;
    case enc_float :
      for (i=0; i<n; i++)
	v(i) = big_endian(read_float(file));
      break;
    case enc_double :
      for (i=0; i<n; i++)
	v(i) = big_endian(read_double(file));
      break;
    default :
      it_warning("SND_In_File::read(): Unsupported encoding!");
      return false;
    }
    return file.good();
  }

  bool SND_In_File::read(vec &v, int n)
  {
    if (!good())
      return false;

    int i;
    
    v.set_size(n, false);
    switch (header.encoding) {
    case enc_linear8 :
      for (i=0; i<n; i++)
	v(i) = read_char(file) / 128.0;
      break;
    case enc_linear16 :
      for (i=0; i<n; i++)
	v(i) = big_endian(read_short(file)) / 32768.0;
      break;
    case enc_float :
      for (i=0; i<n; i++)
	v(i) = big_endian(read_float(file));
      break;
    case enc_double :
      for (i=0; i<n; i++)
	v(i) = big_endian(read_double(file));
      break;
    default :
      it_warning("SND_In_File::read(): Unsupported encoding!");
      return false;
    }
    return file.good();
  }

  //////////////////////////////////////////////////
  //
  // SND_Out_File
  //
  //////////////////////////////////////////////////
  SND_Out_File::SND_Out_File()
  {
  }

  SND_Out_File::SND_Out_File(const char *fname, int rate, data_encoding e)
  {
    open(fname, rate, e);
  }

  bool SND_Out_File::open(const char *fname, int rate, data_encoding e)
  {
    if (file.is_open())
      close();
    file.clear();
    is_valid = false;
    file.open(fname, ios::out | ios::trunc | ios::binary);
    if (!file)
      return false;
    
    header.data_size = 0;
    header.encoding = (unsigned)e;
    header.sample_rate = rate;
    header.channels = 1;
    
    if (!write_header(file))
      return false;

    is_valid = true;
    return true;
  }

  void SND_Out_File::close()
  {
    file.seekp(0, ios::end);
    header.data_size = static_cast<int>(file.tellp()) - sizeof(header);
    write_header(file);
    file.close();
    is_valid = false;
  }
    
  bool SND_Out_File::seek_write(int pos)
  {
    if (!good())
      return false;

    if (pos < 0)
      file.seekp(0, ios::end);
    else
      file.seekp(sizeof(header) + header.channels * sample_size() * pos);
    return true;
  }

  int SND_Out_File::tell_write()
  {
    if (!good())
      return -1;

    return ( static_cast<int>(file.tellp()) - sizeof(header) ) / ( header.channels * sample_size() );
  }

  bool SND_Out_File::write(const vec &v)
  {
    if (!good())
      return false;

    int i;
    
    switch (header.encoding) {
    case enc_linear8 :
      for (i=0; i<v.size(); i++)
	write_char(file, double_to_char(v(i) * 128.0));
      break;
    case enc_linear16 :
      for (i=0; i<v.size(); i++)
	write_short(file, big_endian(double_to_short(v(i) * 32768.0)));
      break;
    case enc_float :
      for (i=0; i<v.size(); i++)
	write_float(file, big_endian((float)v(i)));
      break;
    case enc_double :
      for (i=0; i<v.size(); i++)
	write_double(file, big_endian((double)v(i)));
      break;
    default :
      it_warning("SND_Out_File::write(): Unsupported encoding!");
      return false;
    }

    return file.good();
  }

  //////////////////////////////////////////////////
  //
  // SND_IO_File
  //
  //////////////////////////////////////////////////
  bool SND_IO_File::open(const char *fname)
  {
    if (file.is_open())
      close();
    file.clear();
    is_valid = false;
    file.open(fname, ios::in | ios::out | ios::binary);
    if (!file)
      return false;
    
    if (!read_header(file)) {
      file.close();
      return false;
    }

    if (!seek_read(0) || !seek_write(0)) {
      file.close();
      return false;
    }
    
    is_valid = true;
    return true;
  }

  void SND_IO_File::close()
  {
    write_header(file);
    file.close();
    is_valid = false;
  }
    
  bool snd_read(const char *fname, vec &v)
  {
    SND_In_File file;

    if (!file.open(fname))
      return false;

    return file.read(v);
  }

  bool snd_read(const char *fname, vec &v, int beg, int len)
  {
    SND_In_File file;

    if (!file.open(fname))
      return false;

    file.seek_read(beg);
    return file.read(v, len);
  }

  bool snd_write(const char *fname, const vec &v, int rate, SND_Format::data_encoding e)
  {
    SND_Out_File file;

    if (!file.open(fname, rate, e))
      return false;

    return file.write(v);
  }

} // namespace itpp
