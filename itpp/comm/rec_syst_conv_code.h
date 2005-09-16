/*!
 * \file 
 * \brief Definitions of a Recursive Systematic Convolutional codec class
 * \author Pal Frenger
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
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
 * -------------------------------------------------------------------------
 */

#ifndef REC_SYST_CONV_CODE_H
#define REC_SYST_CONV_CODE_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/binary.h>
#include <itpp/base/specmat.h>
#include <itpp/base/matfunc.h>
#include <itpp/comm/convcode.h>

namespace itpp {

  /*! 
    \ingroup fec
    \brief A Recursive Systematic Convolutional Encoder/Decoder class
  
    The main purpose of this class is its use un the Turbo_Codec class which uses two recursive systematic
    convolutional encoders. It can however be used as a stand alone class. The map_decode
    member function implementation follows the paper "A Turbo Code Tutorial" by William 
    E. Ryan, New Mexico State University. This paper was found on the web and is probably
    unpublished. 
  */
  class Rec_Syst_Conv_Code {
  public:

    //! Class constructor
    Rec_Syst_Conv_Code(void) {}

    //! Class constructor
    virtual ~Rec_Syst_Conv_Code(void) {}

    /*! 
      \brief Set generator polynomials. 

      The generator polynomials are given in Proakis integer form. First generator (gen(0)) is the recursive polynomial.

      \param gen A vector containing the generator polynomials of the RSC Code.
      \param constraint_length The Constraing length of the encoder.
    */
    void set_generator_polynomials(const ivec &gen, int constraint_length);

    /*!
      \brief Sets the channel parameters needed for MAP-decoding.

      \param Ec is the energy per channel symbol
      \param N0 is the single-sided power spectral density of the AWGN on the channel.
    */
    void set_awgn_channel_parameters(double Ec, double N0);

    /*! 
      \brief Set scaling factor for the decoder

      \param in_Lc is the channel reliability factor (i.e. \a Lc = 4 x sqrt( \a Ec ) / \a N0)
    */
    void set_scaling_factor(double in_Lc);

    /*!
      \brief Encode a binary vector of inputs and also adds a tail of \a K-1 zeros to force the encoder into the zero state.

      The encoder remembers that the trellis is terminated in the zero state at the end of the input block. This is then
      utilized by the decoder when going through the trellis in the reverse direction. The tailbits used are returned in \a tail.
      Parity bits for both the input part and the tail part of the data are returned in the matrix \a parity bits.
    */
    void encode_tail(const bvec &input, bvec &tail, bmat &parity_bits);

    /*!
      \brief Encode a binary vector of inputs starting from zero state without adding of a tail
    */
    void encode(const bvec &input, bmat &parity_bits);

    /*!
      \brief Maximum Aposteriori Probability symbol-by-symbol decoder.

      The extrinsic_input is the a priori information on each systematic bit. If no a priori information is availiable, this vector should
      contain only zeros. The extrinsic_output term may be passed to a subsequent decoder in a Turbo
      scheme. The decision variable is \code L = Lc*rec_systematic + extrinsic_output \endcode where \code Lc = 4*sqrt(Ec)/N0 \endcode 

      \param rec_systematic Including both systematic bits and tail bits (if any)
      \param rec_parity Matrix including all parity bits from all polynomials as well as parity bits from the tail (if terminated)
      \param extrinsic_input For all systematic bits
      \param extrinsic_output For all systematic bits
      \param set_terminated Equal to \a true if the trellis was terminated by the encoder and false otherwise 
    */
    virtual void map_decode(const vec &rec_systematic, const mat &rec_parity, const vec &extrinsic_input, vec &extrinsic_output, 
			    bool set_terminated = false);

    /*!
      \brief Log domain implementation of the Maximum Aposteriori Probability symbol-by-symbol decoder.

      The extrinsic_input is the a priori information on each systematic bit. If no a priori information is availiable, this vector should
      contain only zeros. The extrinsic_output term may be passed to a subsequent decoder in a Turbo
      scheme. The decision variable is \code L = Lc*rec_systematic + extrinsic_output \endcode where \code Lc = 4*sqrt(Ec)/N0 \endcode 

      \param rec_systematic Including both systematic bits and tail bits (if any)
      \param rec_parity Matrix including all parity bits from all polynomials as well as parity bits from the tail (if terminated)
      \param extrinsic_input For all systematic bits
      \param extrinsic_output For all systematic bits
      \param set_terminated Equal to \a true if the trellis was terminated by the encoder and false otherwise 
      \param metric May be "LOGMAP" or "LOGMAX" (default)
    */
    virtual void log_decode(const vec &rec_systematic, const mat &rec_parity, const vec &extrinsic_input, vec &extrinsic_output,
			    bool set_terminated = false, string metric = "LOGMAX");

    /*!
      \brief Special Log-MAP/Log-MAX decoder implementation for \a n = 2

      \param rec_systematic Including both systematic bits and tail bits (if any)
      \param rec_parity Matrix including all parity bits from all polynomials as well as parity bits from the tail (if terminated)
      \param extrinsic_input For all systematic bits
      \param extrinsic_output For all systematic bits
      \param in_terminated Equal to \a true if the trellis was terminated by the encoder and false otherwise 
      \param metric May be "LOGMAP" or "LOGMAX" (default)
    */
    virtual void log_decode_n2(const vec &rec_systematic, const vec &rec_parity, const vec &extrinsic_input, vec &extrinsic_output, 
			       bool in_terminated = false, string metric = "LOGMAX");

  private:

    //! Used for precalculations of the trellis state transitions
    int calc_state_transition(const int instate, const int input, ivec &parity);

    int n, K, m;
    ivec gen_pol, gen_pol_rev;
    int encoder_state, Nstates;
    double rate, Lc;
    imat state_trans, output_parity, rev_state_trans, rev_output_parity;
    bool terminated;
    mat gamma, alpha, beta;
    vec denom;
    double ln2;

  };

} // namespace itpp

#endif // #ifndef REC_SYST_CONV_CODE_H
