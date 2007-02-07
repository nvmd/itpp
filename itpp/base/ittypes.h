/*!
 * \file
 * \brief IT++ types definitions
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

#ifndef ITTYPES_H
#define ITTYPES_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#ifndef NO_INT_SIZE_CHECK
#if (SIZEOF_SHORT != 2) || (SIZEOF_UNSIGNED_SHORT != 2) \
  || (SIZEOF_INT != 4) || (SIZEOF_UNSIGNED_INT != 4)
#  error								\
  This platform uses different sizes for "short" and "int" standard	\
  types than expected 2 and 4 bytes, respectively. This causes		\
  incompatibilities of some parts of IT++ with most of 32- and 64-bit	\
  platforms. Especially binary I/O operations will be incompatible.	\
  Please report this problem to IT++ developers. If you are OK with it	\
  you can add "-DNO_INT_SIZE_CHECK" to your CPPFLAGS and recompile the	\
  library.
#endif
#endif // ifndef NO_INT_SIZE_CHECK

#if defined(HAVE_STDINT_H)
#  include <stdint.h>
#elif defined(HAVE_INTTYPES_H)
#  include <inttypes.h>
#else

// Common typedefs for most 32- and 64-bit architechures
typedef signed char             int8_t;	    //!< 8-bit signed integer
typedef unsigned char           uint8_t;    //!< 8-bit unsigned integer
typedef signed short            int16_t;    //!< 16-bit signed integer
typedef unsigned short          uint16_t;   //!< 16-bit unsigned integer
typedef signed int              int32_t;    //!< 32-bit signed integer
typedef unsigned int            uint32_t;   //!< 32-bit unsigned integer

#if defined(_MSC_VER)
typedef __int64                 int64_t;    //!< 64-bit signed integer
typedef unsigned __int64        uint64_t;   //!< 64-bit unsigned integer
#elif (SIZEOF_LONG == 8) && (SIZEOF_UNSIGNED_LONG == 8)
typedef signed long             int64_t;    //!< 64-bit signed integer
typedef unsigned long           uint64_t;   //!< 64-bit unsigned integer
#elif (SIZEOF_LONG_LONG == 8) && (SIZEOF_UNSIGNED_LONG_LONG == 8)
typedef signed long long        int64_t;    //!< 64-bit signed integer
typedef unsigned long long      uint64_t;   //!< 64-bit unsigned integer
#else
#  error						\
  64-bit integer type not detected on this platform.	\
  Please report the problem to IT++ developers.
#endif // defined(_MSC_VER)

#endif // defined(HAVE_STDINT_H)

#endif /* ITTYPES_H */
