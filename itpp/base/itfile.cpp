/*!
 * \file 
 * \brief Implementation of classes for the IT++ file format
 * \author Tony Ottosson and Tobias Ringstrom
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


#include <itpp/base/itfile.h>


namespace itpp {

  char it_file_base::file_magic[4] = { 'I', 'T', '+', '+' };
  char it_file_base::file_version = 2;

  it_ifile::it_ifile()
  {
  }

  it_ifile::it_ifile(const std::string &name)
  {
    open(name);
  }

  void it_ifile::open(const std::string &name)
  {
    if (!exist(name))
      it_error("File does not exist");

    s.open_readonly(name);

    if (!read_check_file_header()) {
      s.close();
     it_error("Corrupt file (Not an it-file)");
    }

  }

  void it_ifile::close()
  {
    s.close();
  }

  bool it_ifile::seek(const std::string &name)
  {
    data_header h;
    std::streampos p;

    s.clear();
    s.seekg(sizeof(file_header));

    while (true) {
      p = s.tellg();
      read_data_header(h);
      if (s.eof()) {
	s.clear();
	return false;
      }
      if (h.type != "" && h.name == name) {
	s.seekg(p);
	break;
      }
      s.seekg(p + static_cast<std::streampos>(h.block_bytes));
    }

    return true;
  }

  bool it_ifile::seek(int n)
  {
    data_header h;
    std::streampos p;

    s.clear();
    s.seekg(sizeof(file_header));
    for (int i=0; i<=n; i++) {
      p = s.tellg(); // changed from tellp() since probably an error
      read_data_header(h);
      if (s.eof()) {
	s.clear();
	return false;
      }
      if (h.type == "")
	i--;
      s.seekg(i==n ? p : p+static_cast<std::streampos>(h.block_bytes));
    }
    return true;
  }

  void it_ifile::info(std::string &name, std::string &type, int &bytes)
  {
    data_header h;
    std::streampos p;

    p = s.tellg(); // changed from tellp()
    read_data_header(h);
    s.seekg(p);
    name = h.name;
    type = h.type;
    bytes = h.data_bytes;
  }

  bool it_ifile::read_check_file_header()
  {
    file_header h;

    memset(&h, 0, sizeof(h)); // Clear the struct
    s.read(reinterpret_cast<char *>(&h), sizeof(h));

    return (memcmp(h.magic, file_magic, 4) == 0 && (h.version <= file_version) );
  }

  void it_ifile::read_data_header(data_header &h)
  {
    std::streampos p=s.tellg();

    s.clear();
    s >> h.endianity;

    if (s.eof())
      return;
    s.set_endianity(static_cast<bfstream_base::endian>(h.endianity));
    s >> h.hdr_bytes;
    s >> h.data_bytes;
    s >> h.block_bytes;
    s >> h.name;
    s >> h.type;
    s.seekg(p + static_cast<std::streampos>(h.hdr_bytes));
  }

  void it_ifile::low_level_read(bin &x)
  {
    s >> x;
  }

  void it_ifile::low_level_read(short &x)
  {
    s >> x;
  }

  void it_ifile::low_level_read(int &x)
  {
    s >> x;
  }

  void it_ifile::low_level_read(float &x)
  {
    s >> x;
  }

  void it_ifile::low_level_read(double &x)
  {
    s >> x;
  }

  void it_ifile::low_level_read(std::complex<float> &x)
  {
    float x_real, x_imag;
    s >> x_real;
    s >> x_imag;
    x = std::complex<float>(x_real, x_imag);
  }

  void it_ifile::low_level_read(std::complex<double> &x)
  {
    double x_real, x_imag;
    s >> x_real;
    s >> x_imag;
    x = std::complex<double>(x_real, x_imag);
  }

  void it_ifile::low_level_read_lo(vec &v)
  {
    int i;
    float val;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++) {
      s >> val;
      v(i) = static_cast<double>(val);
    }
  }

  void it_ifile::low_level_read_hi(vec &v)
  {
    int i;
    double val;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++) {
      s >> val;
      v(i) = static_cast<double>(val);
    }
  }

  void it_ifile::low_level_read(ivec &v)
  {
    int i;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++)
      s >> v(i);
  }

  void it_ifile::low_level_read(bvec &v)
  {
    int i;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++)
      s >> v(i);
  }

  void it_ifile::low_level_read_lo(cvec &v)
  {
    int i;
    float val_real, val_imag;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++) {
      s >> val_real;
      s >> val_imag;
      v(i) = std::complex<double>(val_real, val_imag);
    }
  }

  void it_ifile::low_level_read_hi(cvec &v)
  {
    int i;
    double val_real, val_imag;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++) {
      s >> val_real;
      s >> val_imag;
      v(i) = std::complex<double>(val_real, val_imag);
    }
  }

  void it_ifile::low_level_read(std::string &str)
  {
    int i, j;
    char val;
    str = "";

    s >> i;

    for (j=0; j<i; j++) {
      s >> val;
      str += val;
    }
  }

  void it_ifile::low_level_read_lo(mat &m)
  {
    int i, j;
    float val;

    s >> i >> j;
    m.set_size(i, j, false);
    for (j=0; j<m.cols(); j++)
      for (i=0; i<m.rows(); i++) {
	s >> val;
	m(i,j) = static_cast<double>(val);
      }
  }

  void it_ifile::low_level_read_hi(mat &m)
  {
    int i, j;
    double val;

    s >> i >> j;
    m.set_size(i, j, false);
    for (j=0; j<m.cols(); j++)
      for (i=0; i<m.rows(); i++) {
	s >> val;
	m(i,j) = static_cast<double>(val);
      }
  }

  void it_ifile::low_level_read(imat &m)
  {
    int i, j;

    s >> i >> j;
    m.set_size(i, j, false);
    for (j=0; j<m.cols(); j++)
      for (i=0; i<m.rows(); i++)
	s >> m(i,j);
  }

  void it_ifile::low_level_read(bmat &m)
  {
    int i, j;

    s >> i >> j;
    m.set_size(i, j, false);
    for (j=0; j<m.cols(); j++)
      for (i=0; i<m.rows(); i++)
	s >> m(i,j);
  }

  void it_ifile::low_level_read_lo(cmat &m)
  {
    int i, j;
    float val_real, val_imag;

    s >> i >> j;
    m.set_size(i, j, false);
    for (j=0; j<m.cols(); j++)
      for (i=0; i<m.rows(); i++) {
	s >> val_real;
	s >> val_imag;
	m(i,j) = std::complex<double>(val_real, val_imag);
      }
  }

  void it_ifile::low_level_read_hi(cmat &m)
  {
    int i, j;
    double val_real, val_imag;

    s >> i >> j;
    m.set_size(i, j, false);
    for (j=0; j<m.cols(); j++)
      for (i=0; i<m.rows(); i++) {
	s >> val_real;
	s >> val_imag;
	m(i,j) = std::complex<double>(val_real, val_imag);
      }
  }


  void it_ifile::low_level_read_lo(Array<float> &v)
  {
    int i;
    float val;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++) {
      s >> val;
      v(i) = val;
    }
  }

  void it_ifile::low_level_read_lo(Array<double> &v)
  {
    int i;
    float val;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++) {
      s >> val;
      v(i) = static_cast<double>(val);
    }
  }

  void it_ifile::low_level_read_hi(Array<double> &v)
  {
    int i;
    double val;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++) {
      s >> val;
      v(i) = static_cast<double>(val);
    }
  }

  void it_ifile::low_level_read(Array<int> &v)
  {
    int i;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++)
      s >> v(i);
  }

  void it_ifile::low_level_read(Array<bin> &v)
  {
    int i;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++)
      s >> v(i);
  }

  void it_ifile::low_level_read_lo(Array<std::complex<float> > &v)
  {
    int i;
    float val_real, val_imag;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++) {
      s >> val_real;
      s >> val_imag;
      v(i) = std::complex<float>(val_real, val_imag);
    }
  }

  void it_ifile::low_level_read_lo(Array<std::complex<double> > &v)
  {
    int i;
    float val_real, val_imag;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++) {
      s >> val_real;
      s >> val_imag;
      v(i) = std::complex<double>(val_real, val_imag);
    }
  }

  void it_ifile::low_level_read_hi(Array<std::complex<double> > &v)
  {
    int i;
    double val_real, val_imag;

    s >> i;
    v.set_size(i, false);
    for (i=0; i<v.size(); i++) {
      s >> val_real;
      s >> val_imag;
      v(i) = std::complex<double>(val_real, val_imag);
    }
  }

  it_file::it_file()
  {
    low_prec = false;
    next_name = "";
  }

  it_file::it_file(const std::string &name, bool trunc)
  {
    low_prec = false;
    next_name = "";
    open(name, trunc);
  }

  void it_file::open(const std::string &name, bool trunc)
  {
    if (!exist(name))
      trunc = true;

    s.open(name, trunc);
    it_error_if(!s.is_open(), "Could not open file for writing");

    if (trunc)
      write_file_header();
    else if (!read_check_file_header()) {
      s.close();
      it_error("Corrupt file (Not an it-file)");
    }
  }

  void it_file::close()
  {
    s.close();
  }

  void it_file::flush()
  {
    s.flush();
  }

  void it_file::write_file_header()
  {
    s.write(file_magic, 4);
    s << file_version;
  }

  void it_file::write_data_header(const std::string &type, uint32_t size)
  {
    if (next_name == "")
      it_error("Try to write without a name");
    write_data_header(type, next_name, size);
    next_name = "";
  }

  void it_file::write_data_header(const std::string &type, 
				  const std::string &name, uint32_t size)
  {
    data_header h1, h2;
    std::streampos p;
    int availpos=0;
    bool removed=false;
    int skip;

    h1.endianity = s.get_native_endianity();
    h1.hdr_bytes = 1 + 3*4 + type.size()+1 + name.size()+1;
    h1.data_bytes = size;
    h1.block_bytes = h1.hdr_bytes + h1.data_bytes;
    h1.name = name;
    h1.type = type;

    if (exists(name))
      remove();

    // Try to find an empty space
    s.clear();
    s.seekg(sizeof(file_header));
    while (true) {
      p = s.tellp();
      read_data_header(h2);
      if (s.eof()) {
	s.clear();
	break;
      }
      skip = h2.block_bytes;
      if (h2.type != "" && h2.name == name) {
	s.seekg(p);
	remove();
	s.seekg(p);
	read_data_header(h2);
	removed = true;
	if (availpos != 0)
	  break;
      }
      if (availpos == 0) {
	if (h2.type == "" && h2.block_bytes >= h1.block_bytes) {
	  h1.block_bytes = h2.block_bytes;
	  availpos = p;
	}
	else if (h2.block_bytes-h2.hdr_bytes-h2.data_bytes >= h1.block_bytes) {
	  h1.block_bytes = h2.block_bytes - h2.hdr_bytes - h2.data_bytes;
	  h2.block_bytes = h2.hdr_bytes + h2.data_bytes;
	  s.seekp(p);
	  write_data_header_here(h2);
	  availpos = static_cast<int>(p) + h2.block_bytes;
	  if (removed)
	    break;
	}
      }
      s.seekg(p + static_cast<std::streampos>(skip));
    }
    if (availpos != 0)
      s.seekp(availpos);
    else
      s.seekp(0, std::ios::end);

    write_data_header_here(h1);
  }

  void it_file::write_data_header_here(const data_header &h)
  {
    s.set_endianity(static_cast<bfstream_base::endian>(h.endianity));
    s << h.endianity << h.hdr_bytes << h.data_bytes << h.block_bytes << h.name << h.type;
  }

  void it_file::remove(const std::string &name)
  {
    seek(name);
    remove();
  }

  void it_file::remove()
  {
    data_header h;
    std::streampos p;

    p = s.tellp();
    read_data_header(h);
    h.type = "";
    h.name = "";
    h.hdr_bytes = 1 + 3*4 + 1 + 1;
    h.data_bytes = 0;
    s.seekp(p);
    write_data_header_here(h);
    s.seekp(p + static_cast<std::streampos>(h.block_bytes));
  }

  bool it_file::exists(const std::string &name)
  {
    if (seek(name))
      return true;
    else
      return false;
  }

  void it_file::pack()
  {
    it_warning("pack() is not implemented!");
  }

  void it_file::low_level_write(bin x)
  {
    s << x;
  }

  void it_file::low_level_write(short x)
  {
    s << x;
  }

  void it_file::low_level_write(int x)
  {
    s << x;
  }

  void it_file::low_level_write(float x)
  {
    s << x;
  }

  void it_file::low_level_write(double x)
  {
    s << x;
  }

  void it_file::low_level_write(const std::complex<float> &x)
  {
    s << x.real();
    s << x.imag();
  }

  void it_file::low_level_write(const std::complex<double> &x)
  {
    s << x.real();
    s << x.imag();
  }

  void it_file::low_level_write(const vec &v)
  {
    if (get_low_precision()) {
      s << v.size();
      for (int i=0; i<v.size(); i++)
	s << static_cast<float>(v(i));
    }
    else {
      s << v.size();
      for (int i=0; i<v.size(); i++)
	s << static_cast<double>(v(i));
    }
  }

  void it_file::low_level_write(const ivec &v)
  {
    s << v.size();
    for (int i=0; i<v.size(); i++)
      s << v(i);
  }

  void it_file::low_level_write(const bvec &v)
  {
    s << v.size();
    for (int i=0; i<v.size(); i++)
      s << v(i);
  }

  void it_file::low_level_write(const cvec &v)
  {
    if (get_low_precision()) {
      s << v.size();
      for (int i=0; i<v.size(); i++) {
	s << static_cast<float>(v(i).real());
	s << static_cast<float>(v(i).imag());
      }
    }
    else {
      s << v.size();
      for (int i=0; i<v.size(); i++) {
	s << static_cast<double>(v(i).real());
	s << static_cast<double>(v(i).imag());
      }
    }
  }

  void it_file::low_level_write(const std::string &str)
  {
    s << str.size();

    for (int i=0; i<(int)str.size(); i++)
      s << str[i];
  }

  void it_file::low_level_write(const mat &m)
  {
    int i, j;

    if (get_low_precision()) {
      s << m.rows() << m.cols();
      for (j=0; j<m.cols(); j++)
	for (i=0; i<m.rows(); i++)
	  s << static_cast<float>(m(i,j));
    }
    else {
      s << m.rows() << m.cols();
      for (j=0; j<m.cols(); j++)
	for (i=0; i<m.rows(); i++)
	  s << static_cast<double>(m(i,j));
    }
  }

  void it_file::low_level_write(const imat &m)
  {
    int i, j;

    s << m.rows() << m.cols();
    for (j=0; j<m.cols(); j++)
      for (i=0; i<m.rows(); i++)
	s << m(i,j);
  }

  void it_file::low_level_write(const bmat &m)
  {
    int i, j;

    s << m.rows() << m.cols();
    for (j=0; j<m.cols(); j++)
      for (i=0; i<m.rows(); i++)
	s << m(i,j);
  }

  void it_file::low_level_write(const cmat &m)
  {
    int i, j;

    if (get_low_precision()) {
      s << m.rows() << m.cols();
      for (j=0; j<m.cols(); j++)
	for (i=0; i<m.rows(); i++) {
	  s << static_cast<float>(m(i,j).real());
	  s << static_cast<float>(m(i,j).imag());
	}

    }
    else {
      s << m.rows() << m.cols();
      for (j=0; j<m.cols(); j++)
	for (i=0; i<m.rows(); i++) {
	  s << static_cast<double>(m(i,j).real());
	  s << static_cast<double>(m(i,j).imag());
	}
    }
  }

  void it_file::low_level_write(const Array<float> &v)
  {
    s << v.size();
    for (int i=0; i<v.size(); i++)
      s << v(i);
  }

  void it_file::low_level_write(const Array<double> &v)
  {
    if (get_low_precision()) {
      s << v.size();
      for (int i=0; i<v.size(); i++)
	s << static_cast<float>(v(i));
    }
    else {
      s << v.size();
      for (int i=0; i<v.size(); i++)
	s << static_cast<double>(v(i));
    }
  }

  void it_file::low_level_write(const Array<int> &v)
  {
    s << v.size();
    for (int i=0; i<v.size(); i++)
      s << v(i);
  }

  void it_file::low_level_write(const Array<bin> &v)
  {
    s << v.size();
    for (int i=0; i<v.size(); i++)
      s << v(i);
  }

  void it_file::low_level_write(const Array<std::complex<float> > &v)
  {
    s << v.size();
    for (int i=0; i<v.size(); i++) {
      s << v(i).real();
      s << v(i).imag();
    }
  }

  void it_file::low_level_write(const Array<std::complex<double> > &v)
  {
    if (get_low_precision()) {
      s << v.size();
      for (int i=0; i<v.size(); i++) {
	s << static_cast<float>(v(i).real());
	s << static_cast<float>(v(i).imag());
      }
    }
    else {
      s << v.size();
      for (int i=0; i<v.size(); i++) {
	s << static_cast<double>(v(i).real());
	s << static_cast<double>(v(i).imag());
      }
    }
  }

  it_ifile &operator>>(it_ifile &f, bin &x)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "bin")
      f.low_level_read(x);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, short &x)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "int16")
      f.low_level_read(x);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, int &x)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "int32")
      f.low_level_read(x);
    else if (h.type == "int16") {
      short x16;
      f.low_level_read(x16);
      x = x16;
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, double &x)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "float64")
      f.low_level_read(x);
    else if (h.type == "float32") {
      float f32;
      f.low_level_read(f32);
      x = f32;
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, float &x)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "float32")
      f.low_level_read(x);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, std::complex<float> &x)
  {
    it_file::data_header h;

    f.read_data_header(h);

    if (h.type == "float32_complex") {
      std::complex<float> f32_c;
      f.low_level_read(f32_c);
      x = f32_c;
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, std::complex<double> &x)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "float64_complex")
      f.low_level_read(x);
    else if (h.type == "float32_complex") {
      std::complex<float> f32_c;
      f.low_level_read(f32_c);
      x = f32_c;
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, vec &v)
  {
    it_ifile::data_header h;

    f.read_data_header(h);
    if (h.type == "fvec")
      f.low_level_read_lo(v);
    else if (h.type == "dvec")
      f.low_level_read_hi(v);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, ivec &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "ivec")
      f.low_level_read(v);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, bvec &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "bvec")
      f.low_level_read(v);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, cvec &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "fcvec")
      f.low_level_read_lo(v);
    else if (h.type == "dcvec")
      f.low_level_read_hi(v);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, std::string &str)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "string")
      f.low_level_read(str);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, mat &m)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "fmat")
      f.low_level_read_lo(m);
    else if (h.type == "dmat")
      f.low_level_read_hi(m);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, imat &m)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "imat")
      f.low_level_read(m);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, bmat &m)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "bmat")
      f.low_level_read(m);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, cmat &m)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "fcmat")
      f.low_level_read_lo(m);
    else if (h.type == "dcmat")
      f.low_level_read_hi(m);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<float> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "fArray")
      f.low_level_read_lo(v);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<double> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "fArray")
      f.low_level_read_lo(v);
    else if (h.type == "dArray")
      f.low_level_read_hi(v);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<int> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "iArray")
      f.low_level_read(v);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<bin> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "bArray")
      f.low_level_read(v);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<std::complex<float> > &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "fcArray")
      f.low_level_read_lo(v);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<std::complex<double> > &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "fcArray")
      f.low_level_read_lo(v);
    else if (h.type == "dcArray")
      f.low_level_read_hi(v);
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<vec> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "vecArray") {
      int n;
      f.low_level_read(n);
      v.set_size(n, false);
      for (int i=0; i<n; i++)
	f.low_level_read_hi(v(i));
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<ivec> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "ivecArray") {
      int n;
      f.low_level_read(n);
      v.set_size(n, false);
      for (int i=0; i<n; i++)
	f.low_level_read(v(i));
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<bvec> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "bvecArray") {
      int n;
      f.low_level_read(n);
      v.set_size(n, false);
      for (int i=0; i<n; i++)
	f.low_level_read(v(i));
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<cvec> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "cvecArray") {
      int n;
      f.low_level_read(n);
      v.set_size(n, false);
      for (int i=0; i<n; i++)
	f.low_level_read_hi(v(i));
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<std::string> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "stringArray") {
      int n;
      f.low_level_read(n);
      v.set_size(n, false);
      for (int i=0; i<n; i++)
	f.low_level_read(v(i));
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<mat> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "matArray") {
      int n;
      f.low_level_read(n);
      v.set_size(n, false);
      for (int i=0; i<n; i++)
	f.low_level_read_hi(v(i));
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<imat> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "imatArray") {
      int n;
      f.low_level_read(n);
      v.set_size(n, false);
      for (int i=0; i<n; i++)
	f.low_level_read(v(i));
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<bmat> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "bmatArray") {
      int n;
      f.low_level_read(n);
      v.set_size(n, false);
      for (int i=0; i<n; i++)
	f.low_level_read(v(i));
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_ifile &operator>>(it_ifile &f, Array<cmat> &v)
  {
    it_file::data_header h;

    f.read_data_header(h);
    if (h.type == "cmatArray") {
      int n;
      f.low_level_read(n);
      v.set_size(n, false);
      for (int i=0; i<n; i++)
	f.low_level_read_hi(v(i));
    }
    else
      it_error("Wrong type");

    return f;
  }

  it_file &operator<<(it_file &f, bin x)
  {
    f.write_data_header("bin", sizeof(bin));
    f.low_level_write(x);

    return f;
  }

  it_file &operator<<(it_file &f, short x)
  {
    f.write_data_header("int16", sizeof(short));
    f.low_level_write(x);

    return f;
  }

  it_file &operator<<(it_file &f, int x)
  {
    f.write_data_header("int32", sizeof(int));
    f.low_level_write(x);

    return f;
  }

  it_file &operator<<(it_file &f, float x)
  {
    f.write_data_header("float32", sizeof(float));
    f.low_level_write(x);

    return f;
  }

  it_file &operator<<(it_file &f, double x)
  {
    f.write_data_header("float64", sizeof(double));
    f.low_level_write(x);

    return f;
  }

  it_file &operator<<(it_file &f, std::complex<float> x)
  {
    f.write_data_header("float32_complex", 2*sizeof(float));
    f.low_level_write(x);

    return f;
  }

  it_file &operator<<(it_file &f, std::complex<double> x)
  {
    f.write_data_header("float64_complex", 2*sizeof(double));
    f.low_level_write(x);

    return f;
  }

  it_file &operator<<(it_file &f, const vec &v)
  {
    if (f.get_low_precision())
      f.write_data_header("fvec", sizeof(int) + v.size() * sizeof(float));
    else
      f.write_data_header("dvec", sizeof(int) + v.size() * sizeof(double));
    f.low_level_write(v);

    return f;
  }

  it_file &operator<<(it_file &f, const ivec &v)
  {
    f.write_data_header("ivec", (1 + v.size()) * sizeof(int));
    f.low_level_write(v);

    return f;
  }

  it_file &operator<<(it_file &f, const bvec &v)
  {
    f.write_data_header("bvec", sizeof(int) + v.size() * sizeof(bin) );
    f.low_level_write(v);

    return f;
  }

  it_file &operator<<(it_file &f, const cvec &v)
  {
    if (f.get_low_precision())
      f.write_data_header("fcvec", sizeof(int) + v.size() * 2 * sizeof(float));
    else
      f.write_data_header("dcvec", sizeof(int) + v.size() * 2 * sizeof(double));
    f.low_level_write(v);

    return f;
  }

  it_file &operator<<(it_file &f, const std::string &str)
  {
    f.write_data_header("string", sizeof(int) + str.size() * sizeof(char) );
    f.low_level_write(str);

    return f;
  }

  it_file &operator<<(it_file &f, const mat &m)
  {
    if (f.get_low_precision())
      f.write_data_header("fmat", 2*sizeof(int) + m.rows()*m.cols()*sizeof(float));
    else
      f.write_data_header("dmat", 2*sizeof(int) + m.rows()*m.cols()*sizeof(double));
    f.low_level_write(m);

    return f;
  }

  it_file &operator<<(it_file &f, const imat &m)
  {
    f.write_data_header("imat", (2 + m.rows()*m.cols()) * sizeof(int));
    f.low_level_write(m);

    return f;
  }

  it_file &operator<<(it_file &f, const bmat &m)
  {
    f.write_data_header("bmat", 2*sizeof(int) + m.rows()*m.cols()*sizeof(bin) );
    f.low_level_write(m);

    return f;
  }

  it_file &operator<<(it_file &f, const cmat &m)
  {
    if (f.get_low_precision())
      f.write_data_header("fcmat", 2*sizeof(int) + 2*m.rows()*m.cols()*sizeof(float));
    else
      f.write_data_header("dcmat", 2*sizeof(int) + 2*m.rows()*m.cols()*sizeof(double));
    f.low_level_write(m);

    return f;
  }

  it_file &operator<<(it_file &f, const Array<float> &v)
  {
    f.write_data_header("fArray", sizeof(int) + v.size() * sizeof(float));
    f.low_level_write(v);

    return f;
  }

  it_file &operator<<(it_file &f, const Array<double> &v)
  {
    if (f.get_low_precision())
      f.write_data_header("fArray", sizeof(int) + v.size() * sizeof(float));
    else
      f.write_data_header("dArray", sizeof(int) + v.size() * sizeof(double));
    f.low_level_write(v);

    return f;
  }

  it_file &operator<<(it_file &f, const Array<int> &v)
  {
    f.write_data_header("iArray", (1 + v.size()) * sizeof(int));
    f.low_level_write(v);

    return f;
  }

  it_file &operator<<(it_file &f, const Array<bin> &v)
  {
    f.write_data_header("bArray", sizeof(int) + v.size() * sizeof(bin) );
    f.low_level_write(v);

    return f;
  }

  it_file &operator<<(it_file &f, const Array<std::complex<float> > &v)
  {
    f.write_data_header("fcArray", sizeof(int) + v.size() * 2 * sizeof(float));
    f.low_level_write(v);

    return f;
  }

  it_file &operator<<(it_file &f, const Array<std::complex<double> > &v)
  {
    if (f.get_low_precision())
      f.write_data_header("fcArray", sizeof(int) + v.size() * 2 * sizeof(float));
    else
      f.write_data_header("dcArray", sizeof(int) + v.size() * 2 * sizeof(double));
    f.low_level_write(v);

    return f;
  }

  it_file &operator<<(it_file &f, const Array<vec> &v)
  {
    int i, sum_l=0;

    // calculate total length of Array
    for (i=0; i<v.size(); i++) {
      sum_l += v(i).size();
    }

    // write header
    f.write_data_header("vecArray", sizeof(int)*(1+v.size()) + sum_l * sizeof(double));

    f.low_level_write(v.size());  // the length of the array

    // write one vector at a time (i.e. size and elements)
    for (i=0; i<v.size(); i++)
      f.low_level_write(v(i));

    return f;
  }

  it_file &operator<<(it_file &f, const Array<ivec> &v)
  {
    int i, sum_l=0;

    // calculate total length of Array
    for (i=0; i<v.size(); i++) {
      sum_l += v(i).size();
    }

    // write header
    f.write_data_header("ivecArray", sizeof(int)*(1+v.size()) + sum_l * sizeof(int));

    f.low_level_write(v.size());  // the length of the array

    // write one vector at a time (i.e. size and elements)
    for (i=0; i<v.size(); i++)
      f.low_level_write(v(i));

    return f;
  }

  it_file &operator<<(it_file &f, const Array<bvec> &v)
  {
    int i, sum_l=0;

    // calculate total length of Array
    for (i=0; i<v.size(); i++) {
      sum_l += v(i).size();
    }

    // write header
    f.write_data_header("bvecArray", sizeof(int)*(1+v.size()) + sum_l * sizeof(bin));

    f.low_level_write(v.size());  // the length of the array

    // write one vector at a time (i.e. size and elements)
    for (i=0; i<v.size(); i++)
      f.low_level_write(v(i));

    return f;
  }

  it_file &operator<<(it_file &f, const Array<cvec> &v)
  {
    int i, sum_l=0;

    // calculate total length of Array
    for (i=0; i<v.size(); i++) {
      sum_l += v(i).size();
    }

    // write header
    f.write_data_header("cvecArray", sizeof(int)*(1+v.size()) + sum_l * sizeof(std::complex<double>));

    f.low_level_write(v.size());  // the length of the array

    // write one vector at a time (i.e. size and elements)
    for (i=0; i<v.size(); i++)
      f.low_level_write(v(i));

    return f;
  }

  it_file &operator<<(it_file &f, const Array<std::string> &v)
  {
    int i, sum_l=0;

    // calculate total length of Array
    for (i=0; i<v.size(); i++) {
      sum_l += v(i).size();
    }

    // write header
    f.write_data_header("stringArray", sizeof(int)*(1+v.size()) + sum_l * sizeof(char));

    f.low_level_write(v.size());  // the length of the array

    // write one vector at a time (i.e. size and elements)
    for (i=0; i<v.size(); i++)
      f.low_level_write(v(i));

    return f;
  }

  it_file &operator<<(it_file &f, const Array<mat> &v)
  {
    int i, sum_l=0;

    // calculate total length of Array
    for (i=0; i<v.size(); i++) {
      sum_l += v(i)._datasize();
    }

    // write header
    f.write_data_header("matArray", sizeof(int)*(1+2*v.size()) + sum_l * sizeof(double));

    f.low_level_write(v.size());  // the length of the array

    // write one matrix at a time (i.e. size and elements)
    for (i=0; i<v.size(); i++)
      f.low_level_write(v(i));

    return f;
  }

  it_file &operator<<(it_file &f, const Array<imat> &v)
  {
    int i, sum_l=0;

    // calculate total length of Array
    for (i=0; i<v.size(); i++) {
      sum_l += v(i)._datasize();
    }

    // write header
    f.write_data_header("imatArray", sizeof(int)*(1+2*v.size()) + sum_l * sizeof(int));

    f.low_level_write(v.size());  // the length of the array

    // write one matrix at a time (i.e. size and elements)
    for (i=0; i<v.size(); i++)
      f.low_level_write(v(i));

    return f;
  }

  it_file &operator<<(it_file &f, const Array<bmat> &v)
  {
    int i, sum_l=0;

    // calculate total length of Array
    for (i=0; i<v.size(); i++) {
      sum_l += v(i)._datasize();
    }

    // write header
    f.write_data_header("bmatArray", sizeof(int)*(1+2*v.size()) + sum_l * sizeof(bin));

    f.low_level_write(v.size());  // the length of the array

    // write one matrix at a time (i.e. size and elements)
    for (i=0; i<v.size(); i++)
      f.low_level_write(v(i));

    return f;
  }

  it_file &operator<<(it_file &f, const Array<cmat> &v)
  {
    int i, sum_l=0;

    // calculate total length of Array
    for (i=0; i<v.size(); i++) {
      sum_l += v(i)._datasize();
    }

    // write header
    f.write_data_header("cmatArray", sizeof(int)*(1+2*v.size()) + sum_l * sizeof(std::complex<double>));

    f.low_level_write(v.size());  // the length of the array

    // write one matrix at a time (i.e. size and elements)
    for (i=0; i<v.size(); i++)
      f.low_level_write(v(i));

    return f;
  }

} // namespace itpp
