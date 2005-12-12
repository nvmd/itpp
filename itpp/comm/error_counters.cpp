/*!
 * \file 
 * \brief Implementation of Bit Error Rate Counter (BERC) and 
 * BLock Error Rate Counter (BLERC) classes
 * \author Pal Frenger and Adam Piatyszek
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

#include <itpp/itconfig.h>
#include <itpp/base/binary.h>
#include <itpp/base/matfunc.h>
#include <itpp/comm/error_counters.h>

namespace itpp { 

  //-----------------------------------------------------------
  // The Bit error rate counter class (BERC)
  //-----------------------------------------------------------

  BERC::BERC(long indelay, long inignorefirst, long inignorelast)
  {
    delay       = indelay;
    ignorefirst = inignorefirst;
    ignorelast  = inignorelast;
    errors      = 0;
    corrects    = 0;
  }

  void BERC::count(const bvec &in1, const bvec &in2)
  {
    long countlength = std::min(in1.length(), in2.length()) - std::abs(delay) 
      - ignorefirst - ignorelast;
  
    if (delay >= 0) {
      for (long i = 0; i < countlength; i++) {
	if (in1(i + ignorefirst) == in2(i + ignorefirst + delay)) {
	  corrects++;
	} 
	else {
	  errors++;
	}
      }
    } 
    else {
      for (long i = 0; i < countlength; i++) {
	if (in1(i + ignorefirst - delay) == in2(i + ignorefirst)) {
	  corrects++;
	} 
	else {
	  errors++;
	}
      }
    }
  }

  void BERC::estimate_delay(const bvec &in1, const bvec &in2, long mindelay, 
			    long maxdelay) 
  { 
    long num, start1, start2;
    long min_input_length = std::min(in1.length(), in2.length());
    long bestdelay = mindelay;
    double correlation;
    double bestcorr = 0; 
    for (long i = mindelay; i < maxdelay; i++) {
      // #ifdef _MSC_VER
      // num = min_input_length - abs(i) - ignorefirst - ignorelast;
      // #else
      num = min_input_length - std::abs(i) - ignorefirst - ignorelast;
      // #endif
      start1 = (i < 0) ? -i : 0;
      start2 = (i > 0) ?  i : 0;
      correlation = fabs(sum(to_vec(elem_mult(in1.mid(start1, num), 
					      in2.mid(start2, num)))));
      if (correlation > bestcorr) {
	bestdelay = i;
	bestcorr  = correlation;
      }
    }
    delay = bestdelay;
  }

  void BERC::report()
  {
    std::cout << std::endl
	      << "==================================" << std::endl
	      << "     Bit Error Counter Report     " << std::endl
	      << "==================================" << std::endl
	      << " Ignore First           = " << ignorefirst << std::endl
	      << " Ignore Last            = " << ignorelast << std::endl
	      << " Delay                  = " << delay << std::endl
	      << " Number of counted bits = " << (errors + corrects) << std::endl
	      << " Number of errors       = " << errors << std::endl
	      << "==================================" << std::endl
	      << " Error rate             = " << static_cast<double>(errors) 
      / (errors + corrects) << std::endl
	      << "==================================" << std::endl << std::endl;
  }

  long BERC::count_errors(const bvec &in1, const bvec &in2, long indelay, 
			  long inignorefirst, long inignorelast)
  {
    long countlength = std::min(in1.length(), in2.length()) - std::abs(indelay) 
      - inignorefirst - inignorelast;
    long errors = 0;
  
    if (indelay >= 0) {
      for (long i = 0; i < countlength; i++) {
	if (in1(i + inignorefirst) != in2(i + inignorefirst + indelay)) {
	  errors++;
	}
      }
    } 
    else {
      for (long i = 0; i < countlength; i++) {
	if (in1(i + inignorefirst - indelay) != in2(i + inignorefirst)) {
	  errors++;
	}
      }
    }
    
    return errors;
  }


  //-----------------------------------------------------------
  // The Block error rate counter class (BERC)
  //-----------------------------------------------------------

  BLERC::BLERC(void): setup_done(false), errors(0), corrects(0) {}

  
  BLERC::BLERC(long inblocksize): setup_done(true), blocksize(inblocksize),
				  errors(0), corrects(0) {}

    
  void BLERC::set_blocksize(long inblocksize, bool clear)
  {
    blocksize = inblocksize;
    if (clear) {
      errors = 0;
      corrects = 0;
    }
    setup_done = true;
  }


  void BLERC::count(const bvec &in1, const bvec &in2)
  {
    it_assert(setup_done == true, 
	      "BLERC::count(): Block size has to be setup before counting errors.");
    long min_input_length = std::min(in1.length(), in2.length());
    it_assert(blocksize <= min_input_length, 
	      "BLERC::count(): Block size must not be longer than input vectors."); 
    
    for (int i = 0; i < (min_input_length / blocksize); i++) {
      CORR = true;
      for (int j = 0; j < blocksize; j++) {
	if (in1(i * blocksize + j) != in2(i * blocksize + j)) {
	  CORR = false;
	  break;
	}
      }
      if (CORR) {
	corrects++;
      } else {
	errors++; 
      }
    }
    
  }

} // namespace itpp
