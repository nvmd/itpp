/*!
 * \file
 * \brief Definitions of audio Audio classes and functions
 * \author Tobias Ringstrom and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
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

#ifndef AUDIOFILE_H
#define AUDIOFILE_H

#include <itpp/base/vec.h>
#include <itpp/base/math/misc.h>
#include <fstream>


namespace itpp
{

//! \cond
#define SND_INFO_LEN 8
//! \endcond

/*!
  \addtogroup audio
*/

/*!
  \brief Base class - do not use this one!
  \ingroup audio

  ACTION: ADD DETAILED DOCUMENTATION FOR THIS CLASS!!!!!!!!!!!
*/
class Audio_File
{
public:
  //!Constructor
  Audio_File();
  //!Destructor
  virtual ~Audio_File() { }

  //!Returns true if everything is OK.
  bool good() { return is_valid && file.good(); }

protected:
  //! ACTION: Add documentation for this protected member
  std::fstream file;
  //! ACTION: Add documentation for this protected member
  bool is_valid;
};

/*!
  \brief Base class for SND reading classes (the .au format)
  \ingroup audio

  ACTION: ADD DETAILED DOCUMENTATION FOR THIS CLASS!!!!!!!!!!!
*/
class SND_Format
{
public:
  //! ACTION: ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  enum data_encoding { enc_unknown  =  0,
                       enc_mulaw8   =  1,
                       enc_alaw8    = 27,
                       enc_linear8  =  2,
                       enc_linear16 =  3,
                       enc_linear24 =  4,
                       enc_linear32 =  5,
                       enc_float    =  6,
                       enc_double   =  7
                     };

  //! ACTION: ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  int samples() const { return header.data_size / sample_size(); }
  //! ACTION: ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  data_encoding encoding() const { return (data_encoding)header.encoding; }
  //! ACTION: ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  int rate() const { return header.sample_rate; }
  //! ACTION: ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  void set_rate(int r) { header.sample_rate = r; }
  //! ACTION: ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  int channels() const { return header.channels; }

protected:

  struct {
    //! Magic number
    unsigned magic;
    //! Size of this header
    unsigned hdr_size;
    //! Length of data (optional)
    unsigned data_size;
    //! Data encoding format
    unsigned encoding;
    //! Samples per second
    unsigned sample_rate;
    //! Number of interleaved channels
    unsigned channels;
    //! Info string
    char info[SND_INFO_LEN];
  } header; //!< Definition of the header structure


  //! ACTION: Add documentation for this protected member
  int sample_size() const;
  //! ACTION: Add documentation for this protected member
  bool read_header(std::istream &f);
  //! ACTION: Add documentation for this protected member
  bool write_header(std::ostream &f);
};

/*!
  \brief A class to read SND-files (the .au format)
  \ingroup audio

  ACTION: ADD DETAILED DOCUMENTATION FOR THIS CLASS!!!!!!!!!!!
*/
class SND_In_File : virtual public Audio_File, virtual public SND_Format
{
public:
  //! Default constructor
  SND_In_File();
  //! Open the file {\em fname}.
  SND_In_File(const char *fname);
  //! Destructor
  virtual ~SND_In_File() { close(); }

  //! Open the file {\em fname}.
  virtual bool open(const char *fname);
  //! Close the file.
  virtual void close();

  //! Go to sample number {\em pos}.
  bool seek_read(int pos);
  //! Return the current sample position in the file.
  int tell_read();

  //! Read the whole file into the vector {\em v}.
  virtual bool read(vec &v);
  //! Read {\em n} samples into the vector {\em v}.
  virtual bool read(vec &v, int n);
};

/*!
  \brief A class to write SND-files (the .au format)
  \ingroup audio

  ACTION: ADD DETAILED DOCUMENTATION FOR THIS CLASS!!!!!!!!!!!
*/
class SND_Out_File : virtual public Audio_File, virtual public SND_Format
{
public:
  //! Constructor
  SND_Out_File();
  //! Open the file {\em fname}.
  SND_Out_File(const char *fname, int rate = 8000, data_encoding e = enc_linear16);
  //! Destructor
  virtual ~SND_Out_File() { close(); }

  //! Open the file {\em fname}.
  bool open(const char *fname, int rate = 8000, data_encoding e = enc_linear16);

  // Old definition. Removed since Sun CC gave a warning
  //virtual bool open(const char *fname, int rate=8000, data_encoding e=enc_linear16);

  //! Close the file.
  virtual void close();

  //! Go to sample number {\em pos}.
  bool seek_write(int pos);
  //! Return the current sample position in the file.
  int tell_write();

  //! Write the vector {\em v}.
  virtual bool write(const vec &v);
};

/*!
  \brief This class is capable of doing both input and output.
  \ingroup audio

  ACTION: ADD DETAILED DOCUMENTATION FOR THIS CLASS!!!!!!!!!!!
*/
class SND_IO_File : public SND_In_File, public SND_Out_File
{
public:
  //! Constructor
  SND_IO_File() { }
  //! Open the file {\em fname}.
  SND_IO_File(const char *fname) { open(fname); }
  //! Destructor
  virtual ~SND_IO_File() { close(); }

  //! Open the file {\em fname}.
  virtual bool open(const char *fname);
  //! Close the file.
  virtual void close();
};

/*
   \brief SAP audio file input class
   \ingroup audio

   ADD DETAILED DOCUMENTATION FOR THIS CLASS!!!!!!!!!!!
*/
/*
  class SAP_In_File : virtual public Audio_File {
  public:
  // Constructor
  SAP_In_File();
  // Open the file {\em fname}.
  SAP_In_File(const char *fname);
  // Destructor
  virtual ~SAP_In_File() { close(); }

  // Open the file {\em fname}.
  virtual bool open(const char *fname);
  // Close the file.
  virtual void close();

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  virtual bool seek_read(int pos);
  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  virtual int tell_read();

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  bool read(vec &v);
  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  bool read(vec &v, int n);

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  const char *get_header() { return header; }

  protected:
  char header[SAP_HEADER_SIZE];
  };
*/

/*
  \brief SAP audio file output class
  \ingroup audio

  ADD DETAILED DOCUMENTATION FOR THIS CLASS!!!!!!!!!!!
*/
/*
  class SAP_Out_File : virtual public Audio_File {
  public:
  // Constructor
  SAP_Out_File();
  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  SAP_Out_File(const char *fname, const char *hdr);
  // Destructor
  virtual ~SAP_Out_File() { close(); }

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  bool open(const char *fname, const char *hdr);

  // Old def. Removed since Sun CC gave warning.
  //virtual bool open(const char *fname, const char *hdr);

  // Close the file
  virtual void close();

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  bool seek_write(int pos);
  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  int tell_write();

  // ADD DOCUMENTATION FOR THIS MEMBER!!!!!!!!!!!
  virtual bool write(const vec &v);
  };
*/

/*
  \brief SAP audio file input and output class
  \ingroup audio

  ADD DETAILED DOCUMENTATION FOR THIS CLASS!!!!!!!!!!!
*/
/*
  class SAP_IO_File : public SAP_In_File, public SAP_Out_File {
  public:
  // Constructor
  SAP_IO_File() { }
  // Open the file {\em fname}.
  SAP_IO_File(const char *fname) { open(fname); }
  // Destructor
  virtual ~SAP_IO_File() { close(); }

  // Open the file {\em fname}.
  virtual bool open(const char *fname);
  // Close the file
  virtual void close();
  };
*/

/*! \addtogroup audio */
//!@{

//! Read raw 16-bin little endian audio data
bool raw16le_read(const char *fname, vec &v);
//! Read raw 16-bin little endian audio data
bool raw16le_read(const char *fname, vec &v, int beg, int len);
//! Write raw 16-bin little endian audio data
bool raw16le_write(const char *fname, const vec &v, bool append = false);

//! Read raw 16-bin big endian audio data
bool raw16be_read(const char *fname, vec &v);
//! Read raw 16-bin big endian audio data
bool raw16be_read(const char *fname, vec &v, int beg, int len);
//! Write raw 16-bin big endian audio data
bool raw16be_write(const char *fname, const vec &v, bool append = false);

//! Read SND audio data
bool snd_read(const char *fname, vec &v);
//! Read SND audio data
bool snd_read(const char *fname, vec &v, int beg, int len);
//! Write SND audio data
bool snd_write(const char *fname, const vec &v, int rate = 8000,
               SND_Format::data_encoding e = SND_Format::enc_linear16);
/*
// Read SAP audio data
bool sap_read(const char *fname, vec &v);
// Read SAP audio data
bool sap_read(const char *fname, vec &v, int beg, int len);
// Write SAP audio data
bool sap_write(const char *fname, const vec &v, const char *hdr);
*/

//! Read binary data and optionally switch endianness
template<typename T>
inline T read_endian(std::istream &s, bool switch_endian = false)
{
  T data;
  int bytes = sizeof(T);
  char *c = reinterpret_cast<char *>(&data);
  if (!switch_endian) {
    s.read(c, bytes);
  }
  else {
    for (int i = bytes - 1; i >= 0; i--)
      s.get(c[i]);
  }
  return data;
}

//! Write binary data and optionally switch endianness
template<typename T>
inline void write_endian(std::ostream &s, T data, bool switch_endian = false)
{
  int bytes = sizeof(T);
  char *c = reinterpret_cast<char *>(&data);
  if (!switch_endian) {
    s.write(c, bytes);
  }
  else {
    for (int i = bytes - 1; i >= 0; i--)
      s.put(c[i]);
  }
}

//!@}

} // namespace itpp

#endif // #ifndef AUDIOFILE_H
