/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2005 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Machine dependent formats for binary storage
  \author Tobias Ringström

  $Revision$

  $Date$
*/

#ifndef __machdep_h
#define __machdep_h

#include <itpp/itconfig.h>

namespace itpp {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#define REV_2(from, to)                                   \
              ((char *)(to))[0] = ((char *)(from))[1];    \
              ((char *)(to))[1] = ((char *)(from))[0];

#define REV_4(from, to)                                   \
              ((char *)(to))[0] = ((char *)(from))[3];    \
              ((char *)(to))[1] = ((char *)(from))[2];    \
              ((char *)(to))[2] = ((char *)(from))[1];    \
              ((char *)(to))[3] = ((char *)(from))[0];

#define REV_8(from, to)                                   \
              ((char *)(to))[0] = ((char *)(from))[7];    \
              ((char *)(to))[1] = ((char *)(from))[6];    \
              ((char *)(to))[2] = ((char *)(from))[5];    \
              ((char *)(to))[3] = ((char *)(from))[4];    \
              ((char *)(to))[4] = ((char *)(from))[3];    \
              ((char *)(to))[5] = ((char *)(from))[2];    \
              ((char *)(to))[6] = ((char *)(from))[1];    \
              ((char *)(to))[7] = ((char *)(from))[0];


  //Typedefs for 32 bit architechures (default)
  typedef signed char        it_s8;
  typedef unsigned char      it_u8;
  typedef signed short       it_s16;
  typedef unsigned short     it_u16;
  typedef signed int         it_s32;
  typedef unsigned int       it_u32;
  typedef float              it_f32;
  typedef double             it_f64;

  //One byte data types are independent of endianity:
  inline it_s8 big_endian(it_s8 x)              { return x; }
  inline it_u8 big_endian(it_u8 x)              { return x; }
  inline it_s8 little_endian(it_s8 x)           { return x; }
  inline it_u8 little_endian(it_u8 x)           { return x; }
  inline void  big_endian(it_s8 x, it_s8 &y)    { y = x; }
  inline void  big_endian(it_u8 x, it_u8 &y)    { y = x; }
  inline void  little_endian(it_s8 x, it_s8 &y) { y = x; }
  inline void  little_endian(it_u8 x, it_u8 &y) { y = x; }

#if defined(sparc)
#define __LITTLE_ENDIAN__
#elif defined(i386) || defined(_M_IX86) || defined(__x86_64__) || defined(alpha) || defined(vms)
#define __BIG_ENDIAN__
#endif

  //Determine the endianity (Little or Big):
#if defined (__LITTLE_ENDIAN__)

  //------------------------------------------------------
  // Little Endian
  //------------------------------------------------------

#define IT_ENDIANITY 1234

  inline it_s16 big_endian(it_s16 x) { return x; }
  inline it_u16 big_endian(it_u16 x) { return x; }
  inline it_s32 big_endian(it_s32 x) { return x; }
  inline it_u32 big_endian(it_u32 x) { return x; }
  inline it_f32 big_endian(it_f32 x) { return x; }
  inline it_f64 big_endian(it_f64 x) { return x; }

  inline it_s16 little_endian(it_s16 x) { it_s16 y; REV_2(&x,&y); return y; }
  inline it_u16 little_endian(it_u16 x) { it_u16 y; REV_2(&x,&y); return y; }
  inline it_s32 little_endian(it_s32 x) { it_s32 y; REV_4(&x,&y); return y; }
  inline it_u32 little_endian(it_u32 x) { it_u32 y; REV_4(&x,&y); return y; }
  inline it_f32 little_endian(it_f32 x) { it_f32 y; REV_4(&x,&y); return y; }
  inline it_f64 little_endian(it_f64 x) { it_f64 y; REV_8(&x,&y); return y; }

  inline void big_endian(it_s16 x, it_s16 &y) { y = x; }
  inline void big_endian(it_u16 x, it_u16 &y) { y = x; }
  inline void big_endian(it_s32 x, it_s32 &y) { y = x; }
  inline void big_endian(it_u32 x, it_u32 &y) { y = x; }
  inline void big_endian(it_f32 x, it_f32 &y) { y = x; }
  inline void big_endian(it_f64 x, it_f64 &y) { y = x; }

  inline void little_endian(it_s16 x, it_s16 &y) { REV_2(&x,&y); }
  inline void little_endian(it_u16 x, it_u16 &y) { REV_2(&x,&y); }
  inline void little_endian(it_s32 x, it_s32 &y) { REV_4(&x,&y); }
  inline void little_endian(it_u32 x, it_u32 &y) { REV_4(&x,&y); }
  inline void little_endian(it_f32 x, it_f32 &y) { REV_4(&x,&y); }
  inline void little_endian(it_f64 x, it_f64 &y) { REV_8(&x,&y); }

  //Additions for some 64 bit architechtures
#if defined(sparc)
  typedef signed long        it_s64;
  typedef unsigned long      it_u64;
  inline it_s64 big_endian(it_s64 x) { return x; }
  inline it_u64 big_endian(it_u64 x) { return x; }
  inline it_s64 little_endian(it_s64 x) { it_s64 y; REV_4(&x,&y); return y; }
  inline it_u64 little_endian(it_u64 x) { it_u64 y; REV_4(&x,&y); return y; }
  inline void big_endian(it_s64 x, it_s64 &y) { y = x; }
  inline void big_endian(it_u64 x, it_u64 &y) { y = x; }
  inline void little_endian(it_s64 x, it_s64 &y) { REV_4(&x,&y); }
  inline void little_endian(it_u64 x, it_u64 &y) { REV_4(&x,&y); }
#endif

#elif defined (__BIG_ENDIAN__)

  //------------------------------------------------------
  // Big Endian
  //------------------------------------------------------

#define IT_ENDIANITY 4321

  //32 bit architechures (default)
  inline it_s16 big_endian(it_s16 x) { it_s16 y; REV_2(&x,&y); return y; }
  inline it_u16 big_endian(it_u16 x) { it_u16 y; REV_2(&x,&y); return y; }
  inline it_s32 big_endian(it_s32 x) { it_s32 y; REV_4(&x,&y); return y; }
  inline it_u32 big_endian(it_u32 x) { it_u32 y; REV_4(&x,&y); return y; }
  inline it_f32 big_endian(it_f32 x) { it_f32 y; REV_4(&x,&y); return y; }
  inline it_f64 big_endian(it_f64 x) { it_f64 y; REV_8(&x,&y); return y; }

  inline it_s16 little_endian(it_s16 x) { return x; }
  inline it_u16 little_endian(it_u16 x) { return x; }
  inline it_s32 little_endian(it_s32 x) { return x; }
  inline it_u32 little_endian(it_u32 x) { return x; }
  inline it_f32 little_endian(it_f32 x) { return x; }
  inline it_f64 little_endian(it_f64 x) { return x; }

  inline void big_endian(it_s16 x, it_s16 &y) { REV_2(&x,&y); }
  inline void big_endian(it_u16 x, it_u16 &y) { REV_2(&x,&y); }
  inline void big_endian(it_s32 x, it_s32 &y) { REV_4(&x,&y); }
  inline void big_endian(it_u32 x, it_u32 &y) { REV_4(&x,&y); }
  inline void big_endian(it_f32 x, it_f32 &y) { REV_4(&x,&y); }
  inline void big_endian(it_f64 x, it_f64 &y) { REV_8(&x,&y); }

  inline void little_endian(it_s16 x, it_s16 &y) { y = x; }
  inline void little_endian(it_u16 x, it_u16 &y) { y = x; }
  inline void little_endian(it_s32 x, it_s32 &y) { y = x; }
  inline void little_endian(it_u32 x, it_u32 &y) { y = x; }
  inline void little_endian(it_f32 x, it_f32 &y) { y = x; }
  inline void little_endian(it_f64 x, it_f64 &y) { y = x; }

  //Additions for some 64 bit architechtures
#if defined(alpha) || defined(__x86_64__)
  typedef signed long        it_s64;
  typedef unsigned long      it_u64;
  inline it_s64 big_endian(it_s64 x) { it_s64 y; REV_4(&x,&y); return y; }
  inline it_u64 big_endian(it_u64 x) { it_u64 y; REV_4(&x,&y); return y; }
  inline it_s64 little_endian(it_s64 x) { return x; }
  inline it_u64 little_endian(it_u64 x) { return x; }
  inline void big_endian(it_s64 x, it_s64 &y) { REV_4(&x,&y); }
  inline void big_endian(it_u64 x, it_u64 &y) { REV_4(&x,&y); }
  inline void little_endian(it_s64 x, it_s64 &y) { y = x; }
  inline void little_endian(it_u64 x, it_u64 &y) { y = x; }
#endif

#else

#error "Could not determine endianity!!!"

#endif 

#undef REV_2
#undef REV_4
#undef REV_8

#endif //DOXYGEN_SHOULD_SKIP_THIS

} //namespace itpp

#endif // __machdep_h

