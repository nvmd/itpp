/*!
 * \file 
 * \brief Definition of classes for the IT++ file format
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

#ifndef ITFILE_H
#define ITFILE_H

#include <itpp/base/vec.h>
#include <itpp/base/array.h>
#include <itpp/base/binfile.h>
#include <itpp/base/ittypes.h>


namespace itpp {

  /*!
    \addtogroup itfile
    \brief The IT++ file format
    \author Tony Ottosson and Tobias Ringstrom

    The IT++ file format is a file format that can be used to save variables to
    files.  These files can also be read an written by Matlab using the m-files
    \c load_it.m and \c save_it.m.

    The class it_ifile is used for reading only, and the class it_file
    can be used for both reading and writing.

    The variables must be given a name when they are saved.  The saving is done in
    two steps.  The first step is to supply a name of the variable to be saved.
    This can be done either by calling the function it_file::seek() or by
    using the helper class Name as follows:

    \code
    vec v("1 2 3");
    it_file f("file.it");
    f << Name("v") << v;
    \endcode

    The reading is done in a similar way:

    \code
    vec v;
    it_ifile f("file.it");
    f >> Name("v") >> v;
    \endcode

    \warning Do no use names that begin with an existing type.
  */

  /*!
    \brief Base class for it_ifile and it_file.
    \ingroup itfile
  */
  class it_file_base {
  public:

    //! Data header structure
    struct data_header {
      //! 0=little, 1=big
      char endianity;
      //! size variables
      uint32_t hdr_bytes, data_bytes, block_bytes;
      //! type=="" means deleted
      std::string name, type;
    };

  protected:

    //! File header structure
    struct file_header {
      //! ACTION: Add documentation
      char magic[4];
      //! ACTION: Add documentation
      char version;
    };
    //! ACTION: Add documentation
    static char file_magic[4];
    //! ACTION: Add documentation
    static char file_version;
  };

  /*!
    \brief The IT++ file format reading class.
    \ingroup itfile

  */
  class it_ifile : public it_file_base {
  public:
    //!Constructor.
    it_ifile();
    //!Constructor. Calls open().
    explicit it_ifile(const std::string &name);
    //!Destructor.
    virtual ~it_ifile() { }
    //!Open a file.  The file must exist.
    void open(const std::string &name);
    //!Close a file.
    virtual void close();
    //!Returns pointer to the underlying \c bfstream used
    bfstream &low_level() { return s; }

    //!Reads and checks the file data header. Returns true if the header is valid and false otherwise.
    bool read_check_file_header();
    //!Read the data header and return the result in the variable \c h
    void read_data_header(data_header &h);
    //!Read a binary value at the current file pointer position
    void low_level_read(bin &x);
    //!Read a short value at the current file pointer position
    void low_level_read(short &x);
    //!Read an integer value at the current file pointer position
    void low_level_read(int &x);
    //Read a long long value at the current file pointer position
    //void low_level_read(long_long &x);
    //!Read a float value at the current file pointer position
    void low_level_read(float &x);
    //!Read a double value at the current file pointer position
    void low_level_read(double &x);
    //!Read a float complex value at the current file pointer position
    void low_level_read(std::complex<float> &x);
    //!Read a double complex value at the current file pointer position
    void low_level_read(std::complex<double> &x);
    //!Read a vector of float values at the current file pointer position
    void low_level_read_lo(vec &v);
    //!Read a vector of double values at the current file pointer position
    void low_level_read_hi(vec &v);
    //Read a vector of long long values at the current file pointer position
    //void low_level_read(llvec &v);
    //!Read a vector of integer values at the current file pointer position
    void low_level_read(ivec &v);
    //!Read a vector of binary values at the current file pointer position
    void low_level_read(bvec &v);
    //!Read a vector of float complex values at the current file pointer position
    void low_level_read_lo(cvec &v);
    //!Read a vector of double complex values at the current file pointer position
    void low_level_read_hi(cvec &v);
    //!Read a string at the current file pointer position
    void low_level_read(std::string &str);
    //!Read a matrix of float values at the current file pointer position
    void low_level_read_lo(mat &m);
    //!Read a matrix of double values at the current file pointer position
    void low_level_read_hi(mat &m);
    //Read a matrix of long long values at the current file pointer position
    //void low_level_read(llmat &m);
    //!Read a matrix of integer values at the current file pointer position
    void low_level_read(imat &m);
    //!Read a matrix of binary values at the current file pointer position
    void low_level_read(bmat &m);
    //!Read a matrix of float complex values at the current file pointer position
    void low_level_read_lo(cmat &m);
    //!Read a matrix of double complex values at the current file pointer position
    void low_level_read_hi(cmat &m);

    //!Read an Array of float values at the current file pointer position
    void low_level_read_lo(Array<float> &v);
    //!Read an Array of float values at the current file pointer position
    void low_level_read_lo(Array<double> &v);
    //!Read an Array of double values at the current file pointer position
    void low_level_read_hi(Array<double> &v);
    //Read an Array of long long values at the current file pointer position
    //void low_level_read(Array<long_long> &v);
    //!Read an Array of integer values at the current file pointer position
    void low_level_read(Array<int> &v);
    //!Read an Array of binary values at the current file pointer position
    void low_level_read(Array<bin> &v);
    //!Read an Array of float complex values at the current file pointer position
    void low_level_read_lo(Array<std::complex<float> > &v);
    //!Read an Array of float complex values at the current file pointer position
    void low_level_read_lo(Array<std::complex<double> > &v);
    //!Read an Array of double complex values at the current file pointer position
    void low_level_read_hi(Array<std::complex<double> > &v);

    //! Find the variable \c  name.
    bool seek(const std::string &name);

    //! Find the variable number \c n.
    bool seek(int n);
    //! Get information about the current variable.
    void info(std::string &name, std::string &type, int &bytes);

  protected:
    //! Protected binary file stream
    bfstream s;
  };

  /*!
    \brief The IT++ file format reading and writing class.
    \ingroup itfile

  */
  class it_file : public it_ifile {
  public:
    //! ACTION: Add documentation for this typedef
    typedef it_file& (*it_manip)(it_file&);

    //! Constructor.
    it_file();

    /*!
      \brief Constructor.

      If the file does not exist it will be created.  If \c trunc is true,
      the file will be truncated.
    */
    explicit it_file(const std::string &name, bool trunc=false);

    //! Destructor.
    virtual ~it_file() { }

    /*!
      \brief Open a file for reading and writing.

      If the file does not exist it will be created.  If \c  trunc is true,
      the file will be truncated.
    */
    void open(const std::string &name, bool trunc=false);

    //! Close the file.
    void close();

    //! Flush the data to disk.
    void flush();

    //! Returns pointer to the underlying \c bfstream used
    bfstream &low_level() { return s; }

    //! Set the precision. Low precision means floats, high means doubles.
    void set_low_precision(bool p=true)  { low_prec = p; }

    //! Get the precision.
    bool get_low_precision() { return low_prec; }

    //! Set the name of the next name to be saved. See also the \c Name class.
    void set_next_name(const std::string &n) { next_name=n; }

    //!Write the header for the \c it_file
    void write_file_header();
    //!Write the data header for a variable, specifying the type and size of the data to follow.
    void write_data_header(const std::string &type, uint32_t size);
    //!Write the data header for a variable, specifying the type, name, and size of the data to follow.
    void write_data_header(const std::string &type, const std::string &name,
			   uint32_t size);
    //!Write a binary value at the current file pointer position
    void low_level_write(bin x);
    //!Write a short value at the current file pointer position
    void low_level_write(short x);
    //!Write an integer value at the current file pointer position
    void low_level_write(int x);
    //Write a long long value at the current file pointer position
    //void low_level_write(long_long x);
    //!Write a float value at the current file pointer position
    void low_level_write(float x);
    //!Write a double value at the current file pointer position
    void low_level_write(double x);
    //!Write a float complex value at the current file pointer position
    void low_level_write(const std::complex<float> &x);
    //!Write a double complex value at the current file pointer position
    void low_level_write(const std::complex<double> &x);
    //!Write a vec at the current file pointer position
    void low_level_write(const vec &v);
    //!Write an ivec at the current file pointer position
    void low_level_write(const ivec &v);
    //!Write a bvec at the current file pointer position
    void low_level_write(const bvec &v);
    //Write a llvec at the current file pointer position
    //void low_level_write(const llvec &v);
    //!Write a cvec at the current file pointer position
    void low_level_write(const cvec &v);
    //!Write a string at the current file pointer position
    void low_level_write(const std::string &str);
    //!Write a mat at the current file pointer position
    void low_level_write(const mat &m);
    //!Write a imat at the current file pointer position
    void low_level_write(const imat &m);
    //!Write a bmat at the current file pointer position
    void low_level_write(const bmat &m);
    //Write a llmat at the current file pointer position
    //void low_level_write(const llmat &m);
    //!Write a cmat at the current file pointer position
    void low_level_write(const cmat &m);
    //!Write a float Array at the current file pointer position
    void low_level_write(const Array<float> &v);
    //!Write a double Array at the current file pointer position
    void low_level_write(const Array<double> &v);
    //!Write a integer Array at the current file pointer position
    void low_level_write(const Array<int> &v);
    //!Write a bin Array at the current file pointer position
    void low_level_write(const Array<bin> &v);
    //Write a long long Array at the current file pointer position
    //void low_level_write(const Array<long_long> &v);
    //!Write a float complex Array at the current file pointer position
    void low_level_write(const Array<std::complex<float> > &v);
    //!Write a double complex Array at the current file pointer position
    void low_level_write(const Array<std::complex<double> > &v);

    //!ACTTION: ADD DOCUMENTATION FOR THIS MEMBER !!!!!!!!
    it_file& operator<<(it_manip func) { return (*func)(*this); }

    //! Removes the variable \c name from the file.
    void remove(const std::string &name);
    //! Returns true if the variable \c name exists in the file.
    bool exists(const std::string &name);
    //! Remove slack space from the file.
    void pack();

  protected:
    //! ACTION: Add documenation for this protected member
    void remove();
    //! ACTION: Add documenation for this protected member
    void write_data_header_here(const data_header &h);

    //! ACTION: Add documenation for this protected member
    bool low_prec;
    //! ACTION: Add documenation for this protected member
    std::string next_name;
  };

  /*!
    \brief Flush operator.
    \ingroup itfile

    Flushes the data. Usage:

    \code
    vec v1("1 2 3"), v2;
    it_file f("file.it");
    f << Name("v") << v1 << flush;
    \endcode
  */
  inline it_file& flush(it_file& f)
    {
      f.flush();
      return f;
    }

  /*!
    \brief Automatic naming when saving.
    \ingroup itfile

    An easy way to give a variable a name when saving.  Usage:

    \code
    vec v1("1 2 3"), v2;
    it_file f("file.it");
    f << Name("v") << v1;
    f >> Name("v") >> v2;
    \endcode
  */
  class Name {
  public:
    //! Constructor
    Name(const std::string &n) : name(n) { }
    //! The name string
    const std::string &name;
  };

  /*! \addtogroup itfile */
  //!@{

  //!Finds the variable \c Name in the \c it_ifile. Returns file pointer for reading.
  inline it_ifile &operator>>(it_ifile &f, const Name &s)
    {
      f.seek(s.name);
      return f;
    }

  //!Finds the variable \c Name in the \c it_file. Returns file pointer for writing.
  inline it_file &operator<<(it_file &f, const Name &s)
    {
      f.set_next_name(s.name);
      return f;
    }

  //!Read the binary variable \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, bin &v);

  //!Read the short variable \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, short &v);

  //!Read the integer variable \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, int &v);

  //!Read the float variable \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, float &v);

  //!Read the double variable \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, double &v);

  //!Read the float complex variable \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, std::complex<float> &v);

  //!Read the double complex variable \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, std::complex<double> &v);

  //!Read the vec \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, vec &v);

  //!Read the ivec \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, ivec &v);

  //!Read the bvec \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, bvec &v);

  //!Read the cvec \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, cvec &v);

  //!Read the string \c str from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, std::string &str);

  //!Read the mat \c m from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, mat &m);

  //!Read the imat \c m from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, imat &m);

  //!Read the bmat \c m from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, bmat &m);

  //!Read the cmat \c m from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, cmat &m);

  //!Read the float Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<float> &v);

  //!Read the double Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<double> &v);

  //!Read the integer Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<int> &v);

  //!Read the binary Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<bin> &v);

  //!Read the float complex Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<std::complex<float> > &v);

  //!Read the double complex Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<std::complex<double> > &v);

  //!Read the vec Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<vec> &v);

  //!Read the ivec Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<ivec> &v);

  //!Read the bvec Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<bvec> &v);

  //!Read the cvec Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<cvec> &v);

  //!Read the string Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<std::string> &v);

  //!Read the mat Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<mat> &v);

  //!Read the imat Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<imat> &v);

  //!Read the bmat Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<bmat> &v);

  //!Read the cmat Array \c v from the \c it_ifile pointer
  it_ifile &operator>>(it_ifile &f, Array<cmat> &v);

  //!Write the binary variable \c x to the \c it_file pointer
  it_file &operator<<(it_file &f, bin x);

  //!Write the short variable \c x to the \c it_file pointer
  it_file &operator<<(it_file &f, short x);

  //!Write the integer variable \c x to the \c it_file pointer
  it_file &operator<<(it_file &f, int x);

  //!Write the float variable \c x to the \c it_file pointer
  it_file &operator<<(it_file &f, float x);

  //!Write the double variable \c x to the \c it_file pointer
  it_file &operator<<(it_file &f, double x);

  //!Write the float complex variable \c x to the \c it_file pointer
  it_file &operator<<(it_file &f, std::complex<float> x);

  //!Write the double complex variable \c x to the \c it_file pointer
  it_file &operator<<(it_file &f, std::complex<double> x);

  //!Write the vec \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const vec &v);

  //!Write the ivec \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const ivec &v);

  //!Write the bvec \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const bvec &v);

  //!Write the cvec \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const cvec &v);

  //!Write the string \c str to the \c it_file pointer
  it_file &operator<<(it_file &f, const std::string &str);

  //!Write the mat \c m to the \c it_file pointer
  it_file &operator<<(it_file &f, const mat &m);

  //!Write the imat \c m to the \c it_file pointer
  it_file &operator<<(it_file &f, const imat &m);

  //!Write the bmat \c m to the \c it_file pointer
  it_file &operator<<(it_file &f, const bmat &m);

  //!Write the cmat \c m to the \c it_file pointer
  it_file &operator<<(it_file &f, const cmat &m);

  //!Write the float Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<float> &v);

  //!Write the double Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<double> &v);

  //!Write the int Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<int> &v);

  //!Write the bin Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<bin> &v);

  //!Write the float complex Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<std::complex<float> > &v);

  //!Write the double complex Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<std::complex<double> > &v);

  //!Write the vec Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<vec> &v);

  //!Write the ivec Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<ivec> &v);

  //!Write the bvec Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<bvec> &v);

  //!Write the cvec Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<cvec> &v);

  //!Write the string Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<std::string> &v);

  //!Write the mat Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<mat> &v);

  //!Write the imat Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<imat> &v);

  //!Write the bmat Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<bmat> &v);

  //!Write the cmat Array \c v to the \c it_file pointer
  it_file &operator<<(it_file &f, const Array<cmat> &v);

  //! Save the variable v in the file name.it as the name name.
  template <class T>
    void it_save_var_as(const T &v, const std::string &name)
    {
      it_file f(name + ".it");
      f << Name(name) << v;
      f.close();
    }

  //! Load the variable v from the file name.it as the name name.
  template <class T>
    void it_load_var_as(T &v, const std::string &name)
    {
      it_ifile f(name + ".it");
      f.seek(name);
      f >> v;
      f.close();
    }

  //! A convenient macro.  Calling it_save_var(M) saves M as 'M' in the file 'M.it'.
#define it_save_var(v) it_save_var_as(v,#v)
  //! A convenient macro.  Calling it_load_var(M) loads M as 'M' in the file 'M.it'.
#define it_load_var(v) it_load_var_as(v,#v)

  //!@}

} // namespace itpp

#endif // #ifndef IT_FILE_H

