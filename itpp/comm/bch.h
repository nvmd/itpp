/*!
 * \file 
 * \brief Definition of a BCH encoder/decoder class
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

#ifndef BCH_H
#define BCH_H

#include <itpp/comm/galois.h>
#include <itpp/comm/channel_code.h>


namespace itpp {

  /*!
    \addtogroup fec
  */

  //---------------------- BCH --------------------------------------

  /*! 
    \ingroup fec
    \brief Class for binary, narrow-sense BCH codes.
  
    The notation used is found in S. B. Wicker, "Error control systems for
    digital communication and storage", Appendix E, Prentice-Hall, 1995.

    Example: 
    \code BCH bch(31,21,2,"3 5 5 1") 
    \endcode 
    uses the generator polynomial
    \f$g(x) = x^{10} + x^9 + x^8 + x^6 + x^5 + x^3 + 1\f$, and is capable of
    correcting 2 errors with \a n = 31 and \a k = 21.
  */
  class BCH : public Channel_Code {
  public:
    //! Initialize a (n,k)-code that can correct t errors
    BCH(int in_n, int in_k, int in_t, ivec genpolynom);

    //! Destructor
    virtual ~BCH(){ }

    //! Encode a bvec of indata
    virtual void encode(const bvec &uncoded_bits, bvec &coded_bits);
    //! Encode a bvec of indata
    virtual bvec encode(const bvec &uncoded_bits);

    //! Decode a bvec of coded data.
    virtual void decode(const bvec &coded_bits, bvec &decoded_bits);
    //! Decode a bvec of coded data
    virtual bvec decode(const bvec &coded_bits);

    // Soft-decision decoding is not implemented
    virtual void decode(const vec &received_signal, bvec &output);
    virtual bvec decode(const vec &received_signal);

    //! Get the code rate
    virtual double get_rate() {return double(k)/double(n); }

    //protected:
  private:
    int n, k, t;
    GFX g;
  };

} // namespace itpp

#endif // #ifndef BCH_H
