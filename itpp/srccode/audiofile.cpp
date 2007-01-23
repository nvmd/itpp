/*!
 * \file 
 * \brief Implementation of audio Audio classes and functions
 * \author Tobias Ringstrom and Adam Piatyszek
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
 * This program is free software; you can redistribute it and/or modfy
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

#include <itpp/srccode/audiofile.h>
#include <itpp/base/converters.h>
#include <itpp/base/it_endian.h>
#include <itpp/base/it_types.h>
#include <iostream>


using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::ios;

namespace itpp {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#define SND_MAGIC    0x2e736e64

#endif // #ifndef DOXYGEN_SHOULD_SKIP_THIS

  inline static int16_t double_to_int16_t(double x)
  {
    x *= 32768.0;
    if (x >= 32767.0)
      return 32767;
    else if (x <= -32768.0)
      return -32768;
    else
      return round_i(x);
  }

  inline static int8_t double_to_int8_t(double x)
  {
    x *= 128.0;
    if (x >= 127.0)
      return 127;
    else if (x <= -128.0)
      return -128;
    else
      return round_i(x);
  }

  bool raw16le_read(const char *fname, vec &v)
  {
    ifstream file(fname, ios::in | ios::binary);
    if (!file)
      return false;

    // Check size of a file
    file.seekg(0, ios::end);
    int size = file.tellg();
    file.seekg(0, ios::beg);

    bool switch_endian = check_big_endianness(); // if BIG_ENDIAN than switch
    int n = size / 2; // int16_t vs. byte
    v.set_size(n, false);
    for (int i = 0; i < n; i++)
      v(i) = read_endian<int16_t, ifstream&>(file, switch_endian) / 32768.0;
    
    return true;
  }

  bool raw16le_read(const char *fname, vec &v, int beg, int len)
  {
    it_assert1(len >= 0, "raw16le_read()");
    ifstream file(fname, ios::in | ios::binary);
    if (!file)
      return false;

    bool switch_endian = check_big_endianness(); // if BIG_ENDIAN than switch
    v.set_size(len, false);
    file.seekg(2 * beg);
    for (int i = 0; i < len; i++)
      v(i) = read_endian<int16_t, ifstream&>(file, switch_endian) / 32768.0;
	
    return true;
  }

  bool raw16le_write(const char *fname, const vec &v, bool append)
  {
    ofstream file(fname, (append ? ios::app | ios::ate : ios::out | ios::trunc) | ios::binary);
    if (!file)
      return false;

    bool switch_endian = check_big_endianness(); // if BIG_ENDIAN than switch
    for (int i = 0; i < v.size(); i++)
      write_endian<int16_t, ofstream&>(file, double_to_int16_t(v(i)), switch_endian);
    
    return true;
  }

  bool raw16be_read(const char *fname, vec &v)
  {
    ifstream file(fname, ios::in | ios::binary);
    if (!file)
      return false;

    // Check size of a file
    file.seekg(0, ios::end);
    int size = file.tellg();
    file.seekg(0, ios::beg);

    bool switch_endian = !check_big_endianness(); // if LITTLE_ENDIAN than switch
    int n = size / 2; // int16_t vs. byte
    v.set_size(n, false);
    for (int i = 0; i < n; i++)
      v(i) = read_endian<int16_t, ifstream&>(file, switch_endian) / 32768.0;
    
    return true;
  }

  bool raw16be_read(const char *fname, vec &v, int beg, int len)
  {
    it_assert1(len >= 0, "raw16le_read()");
    ifstream file(fname, ios::in | ios::binary);
    if (!file)
      return false;

    bool switch_endian = !check_big_endianness(); // if LITTLE_ENDIAN than switch
    v.set_size(len, false);
    file.seekg(2 * beg);
    for (int i = 0; i < len; i++)
      v(i) = read_endian<int16_t, ifstream&>(file, switch_endian) / 32768.0;
	
    return true;
  }

  bool raw16be_write(const char *fname, const vec &v, bool append)
  {
    ofstream file(fname, (append ? ios::app | ios::ate : ios::out | ios::trunc) | ios::binary);
    if (!file)
      return false;

    bool switch_endian = !check_big_endianness(); // if LITTLE_ENDIAN than switch
    for (int i = 0; i < v.size(); i++)
      write_endian<int16_t, ofstream&>(file, double_to_int16_t(v(i)), switch_endian);
	
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

  bool SND_Format::read_header(std::istream &f)
  {
    bool switch_endian = !check_big_endianness(); // if LITTLE_ENDIAN than switch
    f.seekg(0);
    header.magic = read_endian<unsigned int, istream&>(f, switch_endian);
    header.hdr_size = read_endian<unsigned int, istream&>(f, switch_endian);
    header.data_size = read_endian<unsigned int, istream&>(f, switch_endian);
    header.encoding = read_endian<unsigned int, istream&>(f, switch_endian);
    header.sample_rate = read_endian<unsigned int, istream&>(f, switch_endian);
    header.channels = read_endian<unsigned int, istream&>(f, switch_endian);
    f.read(header.info, SND_INFO_LEN);
    if (!f || header.magic != SND_MAGIC) {
      std::cerr << header.magic << " != " << SND_MAGIC << std::endl;
      it_warning("SND_Format::read_header(): This is not a .snd file!");
      return false;
    }
    f.seekg(header.hdr_size);

    return f.good();
  }

  bool SND_Format::write_header(std::ostream &f)
  {
    f.seekp(0);
    header.magic = SND_MAGIC;
    header.hdr_size = sizeof(header);
    memset(header.info, 0, SND_INFO_LEN);

    bool switch_endian = !check_big_endianness(); // if LITTLE_ENDIAN than switch
    write_endian<unsigned int, ostream&>(f, header.magic, switch_endian);
    write_endian<unsigned int, ostream&>(f, header.hdr_size, switch_endian);
    write_endian<unsigned int, ostream&>(f, header.data_size, switch_endian);
    write_endian<unsigned int, ostream&>(f, header.encoding, switch_endian);
    write_endian<unsigned int, ostream&>(f, header.sample_rate, switch_endian);
    write_endian<unsigned int, ostream&>(f, header.channels, switch_endian);
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
    
    return ((static_cast<int>(file.tellg()) - sizeof(header)) 
	    / (header.channels * sample_size()));
  }

  bool SND_In_File::read(vec &v)
  {
    if (!good())
      return false;
    
    int i, n;
    
    n = samples();
    v.set_size(n, false);
    seek_read(0);
    
    bool switch_endian = !check_big_endianness(); // if LITTLE_ENDIAN than switch
    switch (header.encoding) {
    case enc_linear8 :
      for (i=0; i<n; i++)
	v(i) = read_endian<int8_t, istream&>(file, switch_endian) / 128.0;
      break;
    case enc_linear16 :
      for (i=0; i<n; i++)
	v(i) = read_endian<int16_t, istream&>(file, switch_endian) / 32768.0;
      break;
    case enc_float :
      for (i=0; i<n; i++)
	v(i) = read_endian<float, istream&>(file, switch_endian);
      break;
    case enc_double :
      for (i=0; i<n; i++)
	v(i) = read_endian<double, istream&>(file, switch_endian);
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
    
    bool switch_endian = !check_big_endianness(); // if LITTLE_ENDIAN than switch
    v.set_size(n, false);
    switch (header.encoding) {
    case enc_linear8 :
      for (i=0; i<n; i++)
	v(i) = read_endian<int8_t, istream&>(file, switch_endian) / 128.0;
      break;
    case enc_linear16 :
      for (i=0; i<n; i++)
	v(i) = read_endian<int16_t, istream&>(file, switch_endian) / 32768.0;
      break;
    case enc_float :
      for (i=0; i<n; i++)
	v(i) = read_endian<float, istream&>(file, switch_endian);
      break;
    case enc_double :
      for (i=0; i<n; i++)
	v(i) = read_endian<double, istream&>(file, switch_endian);
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
    header.encoding = static_cast<unsigned int>(e);
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

    return ((static_cast<int>(file.tellp()) - sizeof(header)) 
	    / (header.channels * sample_size()));
  }

  bool SND_Out_File::write(const vec &v)
  {
    if (!good())
      return false;

    int i;
    
    bool switch_endian = !check_big_endianness(); // if LITTLE_ENDIAN than switch
    switch (header.encoding) {
    case enc_linear8 :
      for (i=0; i<v.size(); i++)
	write_endian<int8_t, ostream&>(file, double_to_int8_t(v(i)), switch_endian);
      break;
    case enc_linear16 :
      for (i=0; i<v.size(); i++)
	write_endian<int16_t, ostream&>(file, double_to_int16_t(v(i)), switch_endian);
      break;
    case enc_float :
      for (i=0; i<v.size(); i++)
	write_endian<float, ostream&>(file, static_cast<float>(v(i)), switch_endian);
      break;
    case enc_double :
      for (i=0; i<v.size(); i++)
	write_endian<double, ostream&>(file, static_cast<double>(v(i)), switch_endian);
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
