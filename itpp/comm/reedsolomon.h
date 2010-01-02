/*!
 * \file
 * \brief Definitions of a Reed-Solomon codec class
 * \author Pal Frenger, Steve Peters and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#ifndef REEDSOLOMON_H
#define REEDSOLOMON_H

#include <itpp/base/vec.h>
#include <itpp/comm/galois.h>
#include <itpp/comm/channel_code.h>


namespace itpp
{

//---------------------- Reed-Solomon --------------------------------------

/*!
  \ingroup fec
  \brief Reed-Solomon Codes.

  Uses the Berlkamp-Massey algorithm for decoding as described in: S. B. Wicker,
  "Error Control Systems for digital communication and storage," Prentice Hall.

  The code is \f$2^m\f$ - ary of length \f$2^m-1\f$ capable of correcting \f$t\f$ errors.
*/
class Reed_Solomon : public Channel_Code
{
public:
  //! Class constructor for the \f$2^m\f$ - ary, \f$t\f$ error correcting RS-code.
  Reed_Solomon(int in_m, int in_t, bool sys = false);
  //! Destructor
  virtual ~Reed_Solomon() { }

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
  virtual double get_rate() const { return static_cast<double>(k) / n; }

  //! Dummy assignment operator - MSVC++ warning C4512
  Reed_Solomon & operator=(const Reed_Solomon &) { return *this; }

protected:
  /*! Internal encoder/decoder parameters
   * @{ */
  int m, t, k, n, q;
  /*! @} */
  //! The generator polynomial of the RS code
  GFX g;
  //! Whether or not the code is systematic
  const bool systematic;
};

} // namespace itpp

#endif // #ifndef REEDSOLOMON_H
