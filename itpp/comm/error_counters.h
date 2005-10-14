/*!
 * \file 
 * \brief Definitions of Bit Error Rate Counter (BERC) and 
 * BLock Error Rate Counter (BLERC) classes
 * \author Pal Frenger
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
 
#ifndef ERROR_COUNTERS_H
#define ERROR_COUNTERS_H

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

} // namespace itpp

#endif // #ifndef ERROR_COUNTERS_H

