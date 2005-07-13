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
  \brief Definition of Class for the Extended Golay Code (24,12,8)
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __egolay_h
#define __egolay_h

#include "itpp/base/vec.h"
#include "itpp/base/mat.h"
#include "itpp/comm/channel_code.h"

namespace itpp {

  /*! 
    \ingroup fec
    \brief Extended Golay code (24,12,8).
    \author Tony Ottosson
  
    The code is given in systematic form with the information bits first, followed by
    the parity check bits. The decoder uses the arithmetic decoding algorithm that is
    for example described in Wicker "Error Control Systems for Digital Communication and
    Storage", Prentice Hall, 1995 (page 143).
  */
  class Extended_Golay : public Channel_Code {
  public:
    //! Constructor
    Extended_Golay();
    //! Destructor
    virtual ~Extended_Golay(){ }

    //! Encoder. Will truncate some bits if not \a length = \c integer * 12
    virtual void encode(const bvec &uncoded_bits, bvec &coded_bits);
    //! Encoder. Will truncate some bits if not \a length = \c integer * 12
    virtual bvec encode(const bvec &uncoded_bits);

    //! Decoder. Will truncate some bits if not \a length = \c integer * 24
    virtual void decode(const bvec &coded_bits, bvec &decoded_bits);
    //! Decoder. Will truncate some bits if not \a length = \c integer * 24
    virtual bvec decode(const bvec &coded_bits);

    // Soft-decision decoding is not implemented
    virtual void decode(const vec &received_signal, bvec &output);
    virtual bvec decode(const vec &received_signal);

    //! Get the code rate
    virtual double get_rate() { return 0.5; };

    //! Gets the generator matrix for the code (also the parity check matrix)
    bmat get_G() { return G; }
  private:
    bmat B,G;
  };

} //namespace itpp

#endif // __egolay_h
