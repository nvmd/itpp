/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2002 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Definitions of Bit error rate counter (BERC) and Block error rate counter (BLERC) classes 
  \author Pål Frenger

  $Revision$

  $Date$
*/
 
#ifndef __error_counters_h
#define __error_counters_h

#include <itpp/base/vec.h>

namespace itpp {

  /*! 
    \brief Bit Error Rate Counter (BERC) Class

    Example:
    \code
    #include "itpp/itcomm.h"

    int main() {
    //Initiate the Bit Error Counter
    BERC berc;

    //Initiate a Binary Symetric Channel with cross-over probability 0.1
    BSC binary_symetric_channel(0.1);
  
    bvec transmitted_bits = randb(100);
    bvec received_bits = binary_symetric_channel(transmitted_bits);

    //Count the number of bit errors
    berc.count(transmitted_bits,received_bits);

    cout << "Estimated bit error probability is " << berc.get_errorrate() << endl;
    } 
    \endcode
  */
  class BERC {
  public:
    /*! 
      \brief Constructor for the berc class.

      <ul>
      <li> \a delay is positive if \a in2 is a delayed replica of \a in1 and negative otherwise. </li>
      <li> \a ignorefirst and \a ignorelast may be used if errors in the begining and/or the end is to be ignored.</li>
      </ul>
    */
    BERC(long indelay = 0, long inignorefirst = 0, long inignorelast = 0);
    //! Cumulative error counter
    void  count(const bvec &in1, const bvec &in2);
    //! Run this member function if the delay between \a in1 and \a in2 is unknown.
    void  estimate_delay(const bvec &in1, const bvec &in2, long mindelay = -100, long maxdelay = 100);
    //! Clears the bit error counter
    void  clear() { errors = 0; corrects = 0; }
    //! Writes an error report
    void  report();
    //! Return the \a delay, assumed or estimated, between \a in1 and \a in2.
    long  get_delay() { return delay; }
    //! Returns the counted number of bit errors
    long  get_errors() { return errors; }
    //! Returns the counted number of corectly received bits
    long  get_corrects() { return corrects; }
    //! Returns the estimated bit error rate.
    double get_errorrate() { return double(errors)/double(corrects+errors); }
    /*!
      \brief static function to allow simple and fast count of bit-errors

      Returns the number of errors between in1 and in2. Typical usage:
      \code
      bvec in1 = randb(100);
      bvec in2 = randb(100);
      long errors = BERC::count_errors(in1, in2);
      \endcode
    */
    static long count_errors(const bvec &in1, const bvec &in2, long indelay=0, long inignorefirst=0, long inignorelast=0);

    //protected:
  private:
    long	delay;
    long	ignorefirst;
    long	ignorelast;
    long	errors;
    long	corrects;
  };

  /*!
    \brief Class for counting block error rates.
  
    Use this class to count block errors in binary vectors.
  */
  class BLERC {
  public:
    //! Class constructor
    BLERC(void);
    //! Set the block size
    void set_blocksize(long inblocksize);
    //! Calculate the number of  block errors between \a in1 and \a in2
    void count(const bvec &in1, const bvec &in2);
    //! Clear the block error counter
    void clear() { errors = 0; corrects = 0; }
    //! Returns the number of block errors
    long get_errors() { return errors; }
    //! Returns the number of correct blocks
    long get_corrects() { return corrects; }
    //! Returns the block error rate
    double get_errorrate() { return double(errors)/double(corrects+errors); }

    //protected:
  private:
    long blocksize;
    long errors;
    long corrects;
    bool CORR;
  };

} //namespace itpp

#endif // __error_counters_h

