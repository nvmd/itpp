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
  \brief Generic Channel Code class
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef _channel_code_h
#define _channel_code_h

#include "itpp/base/vec.h"
#include "itpp/comm/modulator.h"


namespace itpp {

  /*!
    \addtogroup fec
  */

  //---------------------- BCH --------------------------------------

  /*! 
    \ingroup fec
    \brief Generic Channel Code class
  
  */
  class Channel_Code {
  public:
    //! 
    Channel_Code(){}
    //! 
    //virtual ~Channel_Code(){}

    //! Encode a bvec of input
    virtual void encode(const bvec &uncoded_bits, bvec &coded_bits) = 0;
    //! Encode a bvec of input
    virtual bvec encode(const bvec &uncoded_bits) = 0;

    //! Decode a bvec of coded data.
    virtual void decode(const bvec &codedbits, bvec &decoded_bits) = 0;
    //! Decode a bvec of coded data
    //virtual bvec decode(const bvec &coded_bits);
    virtual bvec decode(const bvec &coded_bits) = 0;

    //! Decode a vec of received data.
    virtual void decode(const vec &received_signal, bvec &decoded_bits) = 0;
    //! Decode a vec of received data
    //virtual bvec decode(const vec &received_signal);
    virtual bvec decode(const vec &received_signal) = 0;

    //! Get the code rate
    virtual double get_rate() = 0;
  };


  /*! 
    \ingroup fec
    \brief Dummy Channel Code class

    A dummy code class. Uncoded output.
  */
  class Dummy_Code : public Channel_Code {
  public:
    //! 
      Dummy_Code() {}
    //! 
    virtual ~Dummy_Code(){}

    //! Encode a bvec of input
    virtual void encode(const bvec &uncoded_bits, bvec &coded_bits) { coded_bits = uncoded_bits; }
    //! Encode a bvec of input
    virtual bvec encode(const bvec &uncoded_bits) { return uncoded_bits; }

    //! Decode a bvec of coded data.
    virtual void decode(const bvec &coded_bits, bvec &decoded_bits) { decoded_bits = coded_bits; }
    //! Decode a bvec of coded data.
    virtual bvec decode(const bvec &coded_bits) { return coded_bits; }

    //! Decode a vec of received data. Assumes soft input (BPSK modulated)
    virtual void decode(const vec &received_signal, bvec &decoded_bits) { BPSK bpsk; bpsk.demodulate_bits(received_signal, decoded_bits); }
    //! Decode a vec of received data. Assumes soft input (BPSK modulated)
    virtual bvec decode(const vec &received_signal) { bvec out; decode(received_signal,out); return out; }

    //! Get the code rate
    virtual double get_rate() { return 1.0; }
  };



} //namespace itpp

#endif // _channel_code_h
