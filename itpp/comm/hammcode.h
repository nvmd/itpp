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
  \brief Definitions of a Hamming code class
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __hamming_h
#define __hamming_h

#include "itpp/base/vec.h"
#include "itpp/base/mat.h"
#include "itpp/comm/channel_code.h"

namespace itpp {

  /*! 
    \ingroup fec
    \brief Binary Hamming codes
  */
  class Hamming_Code : public Channel_Code {
  public:
    //! Constructor for \c hamming(n,k). n = pow(2,m)-1 and k = pow(2,m)-m-1.
    Hamming_Code(short m);

    //! Destructor
    virtual ~Hamming_Code(){ }

    //! Hamming encoder. Will truncate some bits if not \a length = \c integer * \a k.
    virtual void encode(const bvec &uncoded_bits, bvec &coded_bits);
    //! Hamming encoder. Will truncate some bits if not \a length = \c integer * \a k.
    virtual bvec encode(const bvec &uncoded_bits);

    //! Hamming decoder. Will truncate some bits if not \a length = \c integer * \a n.
    virtual void decode(const bvec &coded_bits, bvec &decoded_bits);
    //! Hamming decoder. Will truncate some bits if not \a length = \c integer * \a n.
    virtual bvec decode(const bvec &coded_bits);

    // Soft-decision decoding is not implemented
    virtual void decode(const vec &received_signal, bvec &output);
    virtual bvec decode(const vec &received_signal);

    //! Get the code rate
    virtual double get_rate() { return (double)k/(double)n; };

    //! Gets the code length \a n.
    short get_n() { return n; };
    //! Gets the number of information bits per code word, \a k.
    short get_k() { return k; };
    //! Gets the parity check matrix for the code.
    bmat get_H() { return H; };
    //! Gets the generator matrix for the code.
    bmat get_G() { return G; };
  private:
    short n, k; 
    bmat H, G;
    void generate_H(void);
    void generate_G(void);
  };

} //namespace itpp

#endif // __hamming_h
