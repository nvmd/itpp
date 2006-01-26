/*!
 * \file 
 * \brief Binary file formats definitions
 * \author Tony Ottosson and Thomas Eriksson
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
 */

#ifndef BINFILE_H
#define BINFILE_H

#include <itpp/base/binary.h>
#include <fstream>


namespace itpp {

  /*!
    \brief Checks if a filename already exists on the disk
    \ingroup itfile
  */
  bool exist(const std::string &name);

  /*!
    \brief Base class for binary file classes
    \ingroup itfile

    This class serves as a base class for the classes \c bofstream, \c bifstream, and \c bfstream. It
    controls the endianity (i.e. the byte order of multibyte numbers on the disk) of the inhereted classes.
  */
  class bfstream_base {
  public:
    /*!
      \brief Definition of the endian data type

      The Endian defines the order in which multibyte numbers are stored in the file.
      The two orders are called "Little Endian" (\c l_endian ) and "Big Endian" (\c b_endian ).

      "Little Endian" means that the low-order byte of the number is stored at the lowest adress
      (i.e. the little end comes first). "Big Endian" means that the high-order byte of the number is
      stored in memory at the highest address (i.e. the big end comes first)
    */
    enum endian { l_endian, b_endian };

    /*!
      \brief Class Constructor

      \param e Defines the endianity of the class. Possible values are \c l_endian for little endian or \c b_endian for big endian. The default value is \c b_endian.
    */
    bfstream_base(endian e=b_endian);

    /*!
      \brief Returns the endianity of the class (\c l_endian or \c b_endian )
    */
    endian get_endianity() const { return endianity; }

    /*!
      \brief Returns the native endianity for this computer architecture (\c l_endian or \c b_endian )

      Intel processors use "Little Endian" byte ordering while e.g. Motorola processors use "Big Endian" byte ordering.
    */
    endian get_native_endianity() const { return native_endianity; }

    /*!
      \brief Set the endianity for this class
    */
    void set_endianity(endian e) { endianity = e; }

    /*!
      \brief Set the endianity of this class to the native endianity for this computer architecture
    */
    void set_native_endianity() { endianity = native_endianity; }

  protected:
    //! The endianity used by this class
    endian endianity;
    //! The native endianity for this computer architecture
    endian native_endianity;
  };

  /*!
    \brief Binary Outfile Class
    \ingroup itfile
  */
  class bofstream : public bfstream_base, public std::ofstream {
  public:
    /*!
      \brief Class constructor that opens a file and sets the endianity

      \param name The name of the file to open
      \param e Defines the endianity of the class. Possible values are \c l_endian for "Little Endian" or \c b_endian for "Big Endian". The default value is \c b_endian.
    */
    bofstream(const std::string &name, endian e=b_endian);

    //! Class Constructor
    bofstream();

    //! Class Destructor
    ~bofstream() { }

    /*!
      \brief Open a file for writing and set the endianity

      \param name The name of the file to open
      \param e Defines the endianity of the class (default value is \c b_endian )
    */
    void open(const std::string &name, endian e=b_endian);

    //! Writes a \c char variable to the binary output file
    bofstream& operator<<(char a);
    //! Writes a \c bin variable to the binary output file
    bofstream& operator<<(const class bin &a);
    //! Writes an \c int variable to the binary output file
    bofstream& operator<<(int a);
    //! Writes an \c unsigned \c int variable to the binary output file
    bofstream& operator<<(unsigned int a);
    //! Writes a \c short variable to the binary output file
    bofstream& operator<<(short a);
    //! Writes an \c unsigned short variable to the binary output file
    bofstream& operator<<(unsigned short a);
    //! Writes a \c float variable to the binary output file
    bofstream& operator<<(float a);
    //! Writes a \c double variable to the binary output file
    bofstream& operator<<(double a);
    // Writes a \c long \c double variable to the binary output file
    //bofstream& operator<<(long double a);
    //! Writes a \c long \c int variable to the binary output file
    bofstream& operator<<(long int a);
    //! Writes a \c unsigned \c long \c int variable to the binary output file
    bofstream& operator<<(unsigned long int a);
    // Writes a \c long_long variable to the binary output file
    //bofstream& operator<<(long_long a);
    //! ACTION: ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!
    bofstream& operator<<(const char *a);
    //! Writes a \c string variable to the binary output file
    bofstream& operator<<(const std::string &a);
  };

  /*!
    \brief Binary Infile Class
    \ingroup itfile
  */
  class bifstream : public bfstream_base, public std::ifstream {
  public:
    /*!
      \brief Class constructor that opens a file and sets the endianity

      \param name The name of the file to open
      \param e Defines the endianity of the class. Possible values are \c l_endian for "Little Endian" or \c b_endian for "Big Endian". The default value is \c b_endian.
    */
    bifstream(const std::string &name, endian e=b_endian);

    //! Class Constructor
    bifstream();

    //! Class Destructor
    ~bifstream() { }

    /*!
      \brief Open a file for reading and set the endianity

      \param name The name of the file to open
      \param e Defines the endianity of the class (default value is \c b_endian )
    */
    void open(const std::string &name, endian e=b_endian);

    //! Returns the length in bytes of the file
    long length();

    //! Reads a \c char variable from the binary input file
    bifstream& operator>>(char &a);
    //! Reads a \c bin variable from the binary input file
    bifstream& operator>>(class bin &a);
    //! Reads an \c int variable from the binary input file
    bifstream& operator>>(int &a);
    //! Reads an \c unsigned \c int variable from the binary input file
    bifstream& operator>>(unsigned int &a);
    //! Reads a \c short \c int variable from the binary input file
    bifstream& operator>>(short int &a);
    //! Reads an \c unsigned \c short \c int variable from the binary input file
    bifstream& operator>>(unsigned short int &a);
    //! Reads a \c float variable from the binary input file
    bifstream& operator>>(float &a);
    //! Reads a \c double variable from the binary input file
    bifstream& operator>>(double &a);
    // Reads a \c long \c double variable from the binary input file
    //bifstream& operator>>(long double &a);
    //! Reads a \c long \c int variable from the binary input file
    bifstream& operator>>(long int &a);
    //! Reads an \c unsigned \c long \c int variable from the binary input file
    bifstream& operator>>(unsigned long int &a);
    // Reads a \c long_long variable from the binary input file
    //bifstream& operator>>(long_long &a);
    //! ACTION: ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!
    bifstream& operator>>(char *a);
    //! Reads a \c string variable from the binary input file
    bifstream& operator>>(std::string &a);
  };

  /*!
    \brief Binary in/out-file Class
    \ingroup itfile
  */
  class bfstream : public bfstream_base, public std::fstream {
  public:
    /*!
      \brief Class constructor that opens a file and sets the endianity

      \param name The name of the file to open
      \param e Defines the endianity of the class. Possible values are \c l_endian for "Little Endian" or \c b_endian for "Big Endian". The default value is \c b_endian.
    */
    bfstream(const std::string &name, endian e=b_endian);

    //! Class Constructor
    bfstream();

    //! Class Destructor
    ~bfstream() { }

    /*!
      \brief Open a file for reading and writing and set the endianity

      \param name The name of the file to open
      \param trunc ACTION: Add documentation for this parameter (default value is \c false)
      \param e Defines the endianity of the class (default value is \c b_endian )
    */
    void open(const std::string &name, bool trunc=false, endian e=b_endian);

    /*!
      \brief Open a file for reading only and set the endianity

      \param name The name of the file to open
      \param e Defines the endianity of the class (default value is \c b_endian )
    */
    void open_readonly(const std::string &name, endian e=b_endian);

    //! Returns the length in bytes of the file
    long length();

    //! Writes a \c char variable to the binary file
    bfstream& operator<<(char a);
    //! Writes a \c bin variable to the binary file
    bfstream& operator<<(const class bin &a);
    //! Writes an \c int variable to the binary file
    bfstream& operator<<(int a);
    //! Writes an \c unsigned \c int variable to the binary file
    bfstream& operator<<(unsigned int a);
    //! Writes a \c short variable to the binary file
    bfstream& operator<<(short a);
    //! Writes an \c unsigned \c short variable to the binary file
    bfstream& operator<<(unsigned short a);
    //! Writes a \c float variable to the binary file
    bfstream& operator<<(float a);
    //! Writes a \c double variable to the binary file
    bfstream& operator<<(double a);
    // Writes a \c long \c double variable to the binary file
    //bfstream& operator<<(long double a);
    //! Writes a \c long \c int variable to the binary file
    bfstream& operator<<(long int a);
    //! Writes an \c unsigned \c long \c int variable to the binary file
    bfstream& operator<<(unsigned long int a);
    // Writes a \c long_long variable to the binary file
    //bfstream& operator<<(long_long a);
    //! ACTION: ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!
    bfstream& operator<<(const char *a);
    //! Writes a \c string variable to the binary file
    bfstream& operator<<(const std::string &a);

    //! Reads a \c char variable from the binary file
    bfstream& operator>>(char &a);
    //! Reads a \c bin variable from the binary file
    bfstream& operator>>(class bin &a);
    //! Reads an \c int variable from the binary file
    bfstream& operator>>(int &a);
    //! Reads an \c unsigned \c int variable from the binary file
    bfstream& operator>>(unsigned int &a);
    //! Reads a \c short \c int variable from the binary file
    bfstream& operator>>(short int &a);
    //! Reads an \c unsigned \c short \c int variable from the binary file
    bfstream& operator>>(unsigned short int &a);
    //! Reads a \c float variable from the binary file
    bfstream& operator>>(float &a);
    //! Reads a \c double variable from the binary file
    bfstream& operator>>(double &a);
    // Reads a \c long \c double variable from the binary file
    //bfstream& operator>>(long double &a);
    //! Reads a \c long \c int variable from the binary file
    bfstream& operator>>(long int &a);
    //! Reads an \c unsigned \c long \c int variable from the binary file
    bfstream& operator>>(unsigned long int &a);
    // Reads a \c long_long variable from the binary file
    //bfstream& operator>>(long_long &a);
    //! ACTION: ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!
    bfstream& operator>>(char *a);
    //! Reads a \c string variable from the binary file
    bfstream& operator>>(std::string &a);
  };

} //namespace itpp

#endif // #ifndef BINFILE_H













