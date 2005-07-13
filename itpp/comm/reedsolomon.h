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
  \brief Definitions of a Reed-Solomon codec class
  \author Pål Frenger

  $Revision$

  $Date$
*/

#ifndef __reedsolomon_h
#define __reedsolomon_h

#include "itpp/base/vec.h"
#include "itpp/comm/galois.h"
#include "itpp/comm/channel_code.h"

namespace itpp {

  //---------------------- Reed-Solomon --------------------------------------

  /*! 
    \ingroup fec
    \brief Reed-Solomon Codes.
  
    Uses the Berlkamp-Massey algorithm for decoding as described in: S. B. Wicker,
    "Error Control Systems for digital communication and storage," Prentice Hall.

    The code is \f$2^m\f$ - ary of length \f$2^m-1\f$ capable of correcting \f$t\f$ errors.
  */
  class Reed_Solomon : public Channel_Code {
  public:
    //! Class constructor for the \f$2^m\f$ - ary, \f$t\f$ error correcting RS-code.
    Reed_Solomon(int in_m, int in_t);
    //! Destructor
    virtual ~Reed_Solomon(){ }

    //! Encoder function.
    virtual void encode(const bvec &uncoded_bits, bvec &coded_bits);
    //! Encoder function.
    virtual bvec encode(const bvec &uncoded_bits);

    //! Decoder function
    virtual void decode(const bvec &coded_bits, bvec &decoded_bits);
    //! Decoder function
    virtual bvec decode(const bvec &coded_bits);

    // Soft-decision decoding is not implemented
    virtual void decode(const vec &received_signal, bvec &output);
    virtual bvec decode(const vec &received_signal);

    //! Gets the rate of the RS-code.
    virtual double get_rate() { return double(k)/double(n); }

  protected:
    //! Internal encoder/decoder parameters
    int m, t, k, n, q;
    //! The generator polynomial of the RS code
    GFX g;
  };

} //namespace itpp

#endif // __reedsolomon_h
