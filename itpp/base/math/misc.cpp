/*!
 * \file
 * \brief Miscellaneous functions - source file
 * \author Tony Ottosson and Adam Piatyszek
 * 
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2006  (see AUTHORS file for a list of contributors)
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

#include <itpp/base/math/misc.h>


namespace itpp { 

  std::string itpp_version(void)
  {
#ifdef PACKAGE_VERSION
    return std::string(PACKAGE_VERSION);
#else
    return std::string("Warning: Version unknown!");
#endif
  }

  bool check_big_endianness()
  {
    int i = 1;
    char *p = (char *) &i;
    if (p[0] == 1) // Lowest address contains the least significant byte
      return false; // LITTLE_ENDIAN
    else
      return true; // BIG_ENDIAN
  }
  
} //namespace itpp
