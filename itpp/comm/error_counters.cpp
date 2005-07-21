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
  \brief Implementation of Bit error rate counter (BERC) and Block error rate counter (BLERC) classes.
  \author Pål Frenger

  $Revision$

  $Date$
*/

#include <itpp/itconfig.h>
#include <itpp/base/binary.h>
#include <itpp/base/matfunc.h>
#include <itpp/comm/error_counters.h>

using std::cout;
using std::endl;

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
    long countlength = std::min(in1.length(),in2.length()) - std::abs(delay) - ignorefirst - ignorelast;
    long i;
  
    if (delay >= 0) {
      for (i=0; i<countlength; i++) {
	if ( (short)(in1(i+ignorefirst)) == (short)(in2(i+ignorefirst+delay)) ) {
	  corrects += 1;
	} else {
	  errors += 1;
	}
      }
    } else {
      for (i=0; i<countlength; i++) {
	if ( (short)(in1(i+ignorefirst-delay)) == (short)(in2(i+ignorefirst)) ) {
	  corrects += 1;
	} else {
	  errors += 1;
	}
      }
    }	
  }	

  void BERC::estimate_delay(const bvec &in1, const bvec &in2, long mindelay, long maxdelay) 
  { 
    long i, length1 = in1.length(), length2 = in2.length(), num, start1, start2;
    long bestdelay = mindelay;
    double bestcorr = 0, correlation;
    for (i=mindelay; i<maxdelay; i++) {
      //#ifdef _MSC_VER
      //num = std::min(length1,length2) - abs(long(i)) - ignorefirst - ignorelast;
      //#else
      num = std::min(length1,length2) - std::abs(long(i)) - ignorefirst - ignorelast;
      //#endif
      start1  = i<0 ? -i : 0;
      start2  = i>0 ?  i : 0;
      correlation = fabs( sum( to_vec( elem_mult( in1.mid(start1,num), in2.mid(start2,num) ) ) ) );
      if (correlation > bestcorr) {
	bestdelay = i;
	bestcorr  = correlation;
      }
    }
    delay = bestdelay;
  }

  void BERC::report()
  {
    cout << endl;
    cout << "==================================" << endl;
    cout << "     Bit Error Counter Report     " << endl;
    cout << "==================================" << endl;
    cout << " Ignore First           = " << ignorefirst << endl;
    cout << " Ignore Last            = " << ignorelast << endl;
    cout << " Delay                  = " << delay << endl;
    cout << " Number of counted bits = " << (errors + corrects) << endl;
    cout << " Number of errors       = " << errors << endl;
    cout << "==================================" << endl;
    cout << " Error rate             = " << double(errors)/double(errors+corrects) << endl;
    cout << "==================================" << endl;
    cout << endl;
  }

  long BERC::count_errors(const bvec &in1, const bvec &in2, long indelay, long inignorefirst, long inignorelast)
  {
    long countlength = std::min(in1.length(),in2.length()) - std::abs(indelay) - inignorefirst - inignorelast;
    long i, err=0;
  
    if (indelay >= 0) {
      for (i=0; i<countlength; i++) {
	if ( (short)(in1(i+inignorefirst)) != (short)(in2(i+inignorefirst+indelay)) ) {
	  err += 1;
	}
      }
    } else {
      for (i=0; i<countlength; i++) {
	if ( (short)(in1(i+inignorefirst-indelay)) != (short)(in2(i+inignorefirst)) ) {
	  err += 1;
	}
      }
    }
    return err;
  }

  //-----------------------------------------------------------
  // The Block error rate counter class (BERC)
  //-----------------------------------------------------------

  BLERC::BLERC(void)
  {
    errors      = 0;
    corrects    = 0;
  }

  void BLERC::set_blocksize(long inblocksize)
  {
    blocksize = inblocksize;
    errors      = 0;
    corrects    = 0;
  }


  void BLERC::count(const bvec &in1, const bvec &in2)
  {
    long countlength = std::min( in1.length()/blocksize, in2.length()/blocksize);
    long i, j;
  
    for (i=0; i<countlength; i++) {
      CORR = true;
      for (j=0; j<blocksize; j++) {
	if (in1(i*blocksize+j)!=in2(i*blocksize+j)) {
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

} //namespace itpp
