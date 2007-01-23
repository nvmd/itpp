/*!
 * \file
 * \brief IT++ endianness handling functions
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

#ifndef IT_ENDIAN_H
#define IT_ENDIAN_H

#include <fstream>

namespace itpp {

  //!\addtogroup miscfunc
  //!@{

  //! Returns machine endianness: big-endian = true; little-endian = false
  inline bool check_big_endianness()
  {
    int i = 1;
    char* p = (char *) &i;
    if (p[0] == 1) // Lowest address contains the least significant byte
      return false; // LITTLE_ENDIAN
    else
      return true; // BIG_ENDIAN
  }


  //! Read binary data and optionally switch endianness
  template<typename T1, typename T2> inline
  T1 read_endian(T2 s, bool switch_endian = false)
  {
    T1 data;
    unsigned int bytes = sizeof(T1);
    char* c = reinterpret_cast<char *>(&data);
    if (!switch_endian) {
      s.read(c, bytes);
    }
    else {
      for (unsigned int i = bytes-1; i >= 0; i--)
	s.get(c[i]);
    }
    return data;
  }
 
  //! Write binary data and optionally switch endianness
  template<typename T1, typename T2> inline
  void write_endian(T2 s, T1 data, bool switch_endian = false)
  {
    unsigned int bytes = sizeof(T1);
    char* c = reinterpret_cast<char *>(&data);
    if (!switch_endian) {
      s.write(c, bytes);
    }
    else {
      for (unsigned int i = bytes-1; i >= 0; i--)
	s.put(c[i]);
    }
  }

  //!@}
  
} // namespace itpp

#endif /* IT_ENDIAN_H */
