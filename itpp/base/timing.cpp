/*!
 * \file
 * \brief Implementation of Timing classes
 * \author Thomas Eriksson, Tony Ottosson and Tobias Ringstrom
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

#include <ctime>
#include <iostream>
#include <cmath>
#ifndef _MSC_VER
#include <sys/time.h>
#endif

#include <itpp/base/timing.h>

namespace itpp { 

  //! Global object for tic and toc functions
  Real_Timer __tic_toc_timer; 

  //----------------------------------------------------------------------------
  //	class Timer
  //----------------------------------------------------------------------------
  Timer::Timer()
  {
    reset();
  }

  void Timer::start(void)
  {
    if (!running) {
      start_time = get_current_time();
      running = true;
    }
  }

  double Timer::stop(void)
  {
    if (running) {
      stop_time = get_current_time();
      elapsed_time += stop_time-start_time;
      running = false;
    }

    return elapsed_time;
  }

  void Timer::reset(double t)
  {
    elapsed_time = t;
    start_time = 0;
    stop_time = 0;
    running = false;
  }

  double Timer::get_time() const
  {
    return running ?
      elapsed_time + get_current_time() - start_time :
      elapsed_time;
  }

  void Timer::tic(void)
  {
    reset();
    start();
  }

  double Timer::toc(void)
  {
    return get_time() ;
  }

  void Timer::toc_print(void)
  {
    std::cout << "Elapsed time = " << get_time() << " seconds" << std::endl;
  }

  //----------------------------------------------------------------------------
  //	class CPU_Timer
  //----------------------------------------------------------------------------
  double CPU_Timer::get_current_time() const
  {
    return static_cast<double>(clock()) / CLOCKS_PER_SEC;
  }

#ifdef _MSC_VER
#include <windows.h>

  /*
    struct timeval {
    time_t         tv_sec;
    long           tv_usec; 
    };
    typedef struct _FILETIME {
    unsigned long dwLowDateTime;
    unsigned long dwHighDateTime;
    } FILETIME;
    void __stdcall GetSystemTimeAsFileTime(FILETIME*);
  */

  int gettimeofday(struct timeval* p, void* tz)
  {
    union {
      long long ns100; /* time since 1 Jan 1601 in 100ns units */
      FILETIME ft;
    } _now;

    GetSystemTimeAsFileTime( &(_now.ft) );
    p->tv_usec=(long)((_now.ns100 / 10LL) % 1000000LL );
    p->tv_sec= (long)((_now.ns100-(116444736000000000LL))/10000000LL); /* time since 1 Jan 1970 */
    return 0;
  }
#endif
  //----------------------------------------------------------------------------
  //	class Real_Timer
  //----------------------------------------------------------------------------
  double Real_Timer::get_current_time() const
  {
#ifdef MINGW
    // gettimeofday() is not defined in sys/time.h when compiling with MinGW.
    // time() only gives 1-sec accuracy instead of 1-microsec accuracy.
    return time(0);
#else
    struct timeval t;
    gettimeofday(&t, 0);
    return t.tv_sec + t.tv_usec * 1.0e-6;
#endif
  }


  void tic()
  {
    __tic_toc_timer.tic();
  }

  double toc()
  {
    return __tic_toc_timer.toc();
  }

  void toc_print()
  {
    __tic_toc_timer.toc_print();
  }

  void pause(double t)
  {
    if (t==-1) {
      std::cout << "(Press enter to continue)" << std::endl;
      getchar();
    } else {
      Real_Timer	T;
      T.start();
      while (T.get_time()<t);
    }
  }

} // namespace itpp
