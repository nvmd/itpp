/*---------------------------------------------------------------------------*
*                                   IT++			             *
*---------------------------------------------------------------------------*
* Copyright (c) 1995-2004 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
\brief Implementation of Timing classes.
\author Thomas Eriksson, Tony Ottosson, and Tobias Ringström

$Revision$

$Date$
*/

#include <ctime>
#include <iostream>
#include <math.h>
#ifndef _MSC_VER
#include <sys/time.h>
#endif

#include <itpp/base/timing.h>

using std::cout;
using std::endl;

namespace itpp { 

    //! Global object for tic and toc functions
    Real_Timer __tic_toc_timer; 

    //--------------------------------------------------------------------------------
    //	class Timer
    //--------------------------------------------------------------------------------
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
            elapsed_time + get_current_time()-start_time :
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
        cout << "Elapsed time = " << get_time() << " seconds" << endl;
    }

    //--------------------------------------------------------------------------------
    //	class CPU_Timer
    //--------------------------------------------------------------------------------
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
            long long ns100; /*time since 1 Jan 1601 in 100ns units */
            FILETIME ft;
        } _now;

        GetSystemTimeAsFileTime( &(_now.ft) );
        p->tv_usec=(long)((_now.ns100 / 10LL) % 1000000LL );
        p->tv_sec= (long)((_now.ns100-(116444736000000000LL))/10000000LL); /* time since 1 Jan 1970 */
        return 0;
    }
#endif
        //--------------------------------------------------------------------------------
        //	class Real_Timer
        //--------------------------------------------------------------------------------
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
            cout << "(Press enter to continue)" << endl ;
            getchar();
        } else {
            Real_Timer	T;
            T.start();
            while (T.get_time()<t);
        }
    }
} //namespace itpp

