/*!
 * \file 
 * \brief Definition of modulator classes
 * \author Tony Ottosson and Adam Piatyszek
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

#ifndef MODULATOR_H
#define MODULATOR_H

#include <iostream>
#include <itpp/base/vec.h>
#include <itpp/base/specmat.h>

namespace itpp {

  /*! \addtogroup modulators
   */

  /*! 
    \ingroup modulators
    \brief Abstract base class for modulators
  */
  class Modulator {
  public:
    //! Constructor
    Modulator() {}
    //! Destructor
    virtual ~Modulator() {}
    //! Returns number of bits per symbol for the modulator. Can be noninteger.
    virtual double bits_per_symbol() const = 0;
    //! Modulation of bits
    virtual void modulate_bits(const bvec &bits, cvec &out) const = 0;
    //! Modulation of bits
    virtual cvec modulate_bits(const bvec &bits) const = 0;
    //! Demodulation of bits
    virtual void demodulate_bits(const cvec &signal, bvec &out) const = 0;
    //! Demodulation of bits
    virtual bvec demodulate_bits(const cvec &signal) const = 0;
    //! Soft demodulator for AWGN channel
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits) const = 0;
    //! Soft demodulator for a known channel in AWGN
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const = 0;
    //! Approximate soft demodulator for AWGN channel
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits) const = 0;
    //! Approximate soft demodulator for a known channel in AWGN
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const = 0;
  };

  
  /*! 
    \ingroup modulators
    \brief A general modulator class for 1-dimensional signal consellations.
  */
  class Modulator_1d {
  public:
    //! Constructor
    Modulator_1d(const vec &insymbols = "0", const ivec &inbitmap = "0");
    //! Destructor
    virtual ~Modulator_1d() { }
    //! Returns number of bits per symbol for the modulator. Can be noninteger.
    virtual double bits_per_symbol() const {return k;}
    //! Modulate function for symbols
    virtual vec modulate(const ivec &symbolnumbers) const;
    //! Demodulate function for symbols
    virtual ivec demodulate(const vec &signal) const;

    //! Modulate function for bits
    virtual vec modulate_bits(const bvec &bits) const;
    //! Demodulate function for bits
    virtual bvec demodulate_bits(const vec &signal) const;

    //! Set the symbol constellation and the corresponding bitmap
    void set(const vec &insymbols, const ivec &inbitmap);

    //! Get the symbol constellation
    vec get_symbols(void) const;
    //! Get the bitmap
    ivec get_bitmap(void) const;
  protected:
    //! Number of bits per modulation symbol
    int k;
    //! Number of modulation symbols
    int M;
    //! Bit mapping table
    ivec bitmap;
    //! Vector of modulation symbols
    vec symbols;
  };

  /*! 
    \ingroup modulators
    \brief General modulator for 2-dimensional signal constellations.

    This class can also perform soft demodulation. To use the soft demodulate member functions the
    received symbols shall equal
    \f[ r_k = c_k \times s_k + n_k, \f]
    where \f$c_k\f$  is the complex channel
    gain, \f$s_k\f$ is the transmitted QAM symbols, and \f$n_k\f$ is the AWGN of the channel (with 
    variance \f$N_0/2\f$ in both the real and the imaginary valued components).
  
    The input samples to the soft demodulate functions shall equal 
    \f[ z_k = \hat{c}_k^{*} \times r_k, \f]
    where \f$\hat{c}_k^{*}\f$ is the conjugate of the channel estimate. This class assumes that the
    channel estimates are perfect when calculating the soft bits.
  
    When these member functions are used together with MAP-based turbo decoding algoritms then the channel
    reliability factor \f$L_c\f$ of the turbo decoder shall be set to 1. The output from these member functions
    can also be used by a Viterbi decoder using an AWGN based metric calculation function. 
  */
  class Modulator_2d {
  public:
    //! Constructor
    Modulator_2d(const cvec &insymbols = zeros_c(1), const ivec &inbitmap = "0");
    //! Destructor
    virtual ~Modulator_2d() {}
    //! Returns number of bits per symbol for the modulator. Can be noninteger.
    virtual double bits_per_symbol() const {return k;}
    //! Modulation of symbols
    virtual cvec modulate(const ivec &symbolnumbers) const;
    //! Demodulation of symbols
    virtual ivec demodulate(const cvec &signal) const;
  
    //! Modulation of bits
    virtual cvec modulate_bits(const bvec &bits) const;
    //! Demodulation of bits
    virtual bvec demodulate_bits(const cvec &signal) const;
  
    /*! 
      \brief Soft demodulator for AWGN channels
    
      This function calculates
      \f[ 
      \log \left( \frac{
      \sum_{s_i \in S_0} \frac{1}{\pi N_0} \exp \left( -\frac{ |z_k - s_i|^2 }{N_0} \right) } {
      \sum_{s_i \in S_1} \frac{1}{\pi N_0} \exp \left( -\frac{ |z_k - s_i|^2 }{N_0} \right) }
      \right)
      \f]
      where \f$s_i \in S_0\f$ denotes a constellation symbol with the i-th bit equal to zero. This function can be used on channels
      where the channel gain \f$c_k = 1\f$.
    
      \param rx_symbols The received noisy constellation symbols
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    void demodulate_soft_bits(const cvec &rx_symbols, double N0, vec &soft_bits);
  
    /*! 
      \brief Soft demodulator for fading channels
    
      This function calculates
      \f[ 
      \log \left( \frac{
      \sum_{s_i \in S_0} \frac{1}{\pi N_0 |c_k|^2} \exp \left( -\frac{ |z_k - |c_k|^2 s_i|^2 }{N_0 |c_k|^2} \right) } {
      \sum_{s_i \in S_1} \frac{1}{\pi N_0 |c_k|^2} \exp \left( -\frac{ |z_k - |c_k|^2 s_i|^2 }{N_0 |c_k|^2} \right) }
      \right)
      \f]
    
      \param rx_symbols The received noisy constellation symbols \f$z_k\f$ (remember that \f$z_k = \hat{c}_k^{*} \times r_k\f$)
      \param chan The complex valued channel values
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    void demodulate_soft_bits(const cvec &rx_symbols, const cvec &chan, const double N0, vec &soft_bits);
  
    /*!
      \brief Approximative soft demodulator for AWGN channels. 
    
      This function is faster and gives allmost no performance degradation compared to the 
      demodulate_soft_bits(const cvec &symbols, vec &soft_bits) function. This function finds for each bit 
      the closest constellation point that have a zero and a one in the corresponding position. Let \f$d_0\f$
      denote the distance to the closest zero point and \f$d_1\f$ denote the distance to the closest one point
      for the corresponding bit respectively. This algorithm then computes
      \f[ \frac{1}{N_0} ( d_1^2 - d_0^2 ) \f]
    
      \param rx_symbols The received noisy constellation symbols
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    void demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits);
  
    /*! 
      \brief Approximative soft demodulator for fading channels. 
    
      This function is faster and gives allmost no performance degradation compared to the 
      demodulate_soft_bits(const cvec &symbols, const cvec &chan, vec &soft_bits) function.
      Let \f$d_0 = | z_k - |c_k|^2 s_0 |\f$  and \f$d_1 = | z_k - |c_k|^2 s_1 |\f$, with 
      \f$s_0\f$ and \f$s_0\f$ denoting the closest constellation points with zero and one in the 
      corresponding bit position respectively. This algorithm then computes
      \f[ \frac{1}{N_0 |c_k|^2} ( d_1^2 - d_0^2 ) \f]
    
      \param rx_symbols The received noisy constellation symbols \f$z_k\f$ (remember that \f$z_k = \hat{c}_k^{*} \times r_k\f$)
      \param chan The complex valued channel values
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    void demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &chan, const double N0, vec &soft_bits);
  
    //! Set the symbol values to use in the modulator.
    void set(const cvec &insymbols, const ivec &inbitmap);
  
    //! Get the symbol values used in the modulator.
    cvec get_symbols() const;
    //! Get the bitmap used in the modulator.
    ivec get_bitmap() const;

  protected:

    //! Number of bits per modulation symbol
    int k;
    //! Number of modulation symbols
    int M;
    //! Bit mapping table
    ivec bitmap;
    //! Vector of modulation symbols
    cvec symbols;
  
    //! This function calculates the soft bit mapping matrices S0 and S1
    void calculate_softbit_matricies(ivec inbitmap);
  
    //! Used by the soft demodulator functions: Matrix where row k contains the constellation points with a zero in bit position k
    imat S0;   
    //! Used by the soft demodulator functions: Matrix where row k contains the constellation points with a one in bit position k
    imat S1;   
    //! Internal protected state variable
    bool soft_bit_mapping_matrices_calculated;
  };

  /*! 
    \ingroup modulators
    \brief BPSK Modulator Class
  
    Symbols used are \f$ \{ 1, -1 \}\f$. Bit mapping: \f$0 \rightarrow 1\f$ and \f$1 \rightarrow -1\f$.
    Example of use:
    \code
    BPSK bpsk;
    bvec bits = "1 0 0 1 1 0 1 0 1 0 1 1 1 0";
    vec symbols = bpsk.modulate_bits(bits);
    \endcode
  */
  class BPSK : virtual public Modulator {
  public:
    //! Constructor
    BPSK() { }
    //! Destructor
    virtual ~BPSK() {}
    //! Returns number of bits per symbol (1 for BPSK)
    virtual double bits_per_symbol() const {return 1;}
    //! Modulate bits into BPSK symbols
    void modulate_bits(const bvec &bits, vec &out) const;
    //! Modulate bits into BPSK symbols and output cvec. Symbols in real part only
    virtual void modulate_bits(const bvec &bits, cvec &out) const;
    //! Modulate bits into BPSK symbols
    virtual cvec modulate_bits(const bvec &bits) const;
    //! Demodulate noisy BPSK symbols into bits
    void demodulate_bits(const vec &signal, bvec &out) const;
    //! Demodulate noisy BPSK symbols into bits
    virtual bvec demodulate_bits(const vec &signal) const;
    /*! 
      \brief Demodulate noisy BPSK symbols into bits.

      Input is cvec but received signals should be multiplied with complex conjugate
      of channel coefficients before demodulation.
    */
    virtual void demodulate_bits(const cvec &signal, bvec &out) const;
    //! Demodulate noisy BPSK symbols into bits
    virtual bvec demodulate_bits(const cvec &signal) const;

    /*! 
      \brief Soft demodulator for AWGN channels
    
      This function calculates the log-MAP estimate assuming equally likely bits transmitted
      \f[
      \log \left( \frac{\Pr(b=0|r)}{\Pr(b=1|r)} \right)
      = \frac{4 r}{N_0}
      \f]
      It is assumed that what is received is \f$r = b + n\f$.
      \param rx_symbols The received noisy constellation symbols, \f$r\f$ (real)
      \param N0 The single sided spectral density of the AWGN noise, \f$n\f$
      \param soft_bits The soft bits calculated using the expression above
    */
    void demodulate_soft_bits(const vec &rx_symbols, const double N0, vec &soft_bits) const;

    /*! 
      \brief Soft demodulator for AWGN channel
    
      This function calculates the log-MAP estimate assuming equally likely bits transmitted
      \f[
      \log \left( \frac{\Pr(b=0|r)}{\Pr(b=1|r)} \right)
      = \frac{4 r}{N_0}
      \f]
      It is assumed that what is received is \f$r = b + n\f$.
      \param rx_symbols The received noisy constellation symbols, \f$r\f$ (complex but symbols in real part)
      \param N0 The single sided spectral density of the AWGN noise, \f$n\f$
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits) const;

    /*! 
      \brief Soft demodulator for a known channel in AWGN      
      This function calculates the log-MAP estimate assuming equally likely bits transmitted
      \f[
      \log \left( \frac{\Pr(b=0|r)}{\Pr(b=1|r)} \right)
      = \frac{ 4 \Re \{r c^{*} \} }{N_0}
      \f]
      It is assumed that what is received is the complex-valued model: \f$r = c b + n\f$.
      \param rx_symbols The received noisy constellation symbols, \f$r\f$ (complex)
      \param channel The channel coefficients, \f$c\f$ (complex)
      \param N0 The single sided spectral density of the AWGN noise, \f$n\f$
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const;

    /*!
      \brief Soft demolulation. Same as the exact given above
     */
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits) const;

    /*!
      \brief Soft demolulation. Same as the exact given above
    */
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const;
  };

  /*! 
    \ingroup modulators
    \brief M-ary PAM modulator
  
    Mary-PAM signals \f$ \{ M-1, \ldots, 3, 1, -1, -3, \ldots, -(M-1) \}\f$.
    Symbol numbering is from right to left in increasing order.
    Gray encoded bitmapping is used.

    The symbols are normalized so that the average energy is 1. That is, normalized with \f$ \sqrt{(M^2-1)/3}\f$.
  */
  class PAM : virtual public Modulator {
  public:
    //! Constructor
    PAM(int Mary) { set_M(Mary); }
    //! Destructor
    virtual ~PAM() { }
    //! Returns number of bits per symbol
    virtual double bits_per_symbol() const {return k;}
    //! Modulate bits into PAM symbols
    void modulate_bits(const bvec &bits, vec &out) const;
    //! Modulate bits into PAM symbols
    virtual void modulate_bits(const bvec &bits, cvec &out) const;
    //! Modulate bits into PAM symbols
    virtual cvec modulate_bits(const bvec &bits) const;
    //! Demodulation of PAM symbols to bits
    void demodulate_bits(const vec &signal, bvec &out) const;
    //! Demodulation of PAM symbols to bits. Symbols in real part
    virtual void demodulate_bits(const cvec &signal, bvec &out) const;
    //! Demodulation of PAM symbols to bits
    virtual bvec demodulate_bits(const cvec &signal) const;

    /*! 
      \brief Soft demodulator for AWGN channels
    
      This function calculates
      \f[ 
      \log \left( \frac{\Pr(b_i=0|r)}{\Pr(b_i=1|r)} \right) =
      \log \left( \frac{
      \sum_{s_i \in S_0} \exp \left( -\frac{ |r - s_i|^2 }{N_0} \right) } {
      \sum_{s_i \in S_1} \exp \left( -\frac{ |r - s_i|^2 }{N_0} \right) }
      \right)
      \f]
      where \f$s_i \in S_0\f$ denotes a constellation symbol with the i-th bit equal to zero.
      This function can be used on channels where the channel gain is \f$c = 1\f$.
    
      \param rx_symbols The received noisy constellation symbols, \f$r\f$ (complex, but symbols are real)
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits) const;

    /*!
      \brief Approximative soft demodulator for AWGN channel
    
      This function is faster and gives almost no performance degradation compared to the 
      demodulate_soft_bits() function. Only include the largest term in nominator and denominator when calculating
      the log-likelihood function.
    */
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits) const;


    /*! 
      \brief Soft demodulator for a known channel in AWGN
    
      This function calculates
      \f[ 
      \log \left( \frac{\Pr(b_i=0|r)}{\Pr(b_i=1|r)} \right) =
      \log \left( \frac{
      \sum_{s_i \in S_0} \exp \left( -\frac{ |r - c s_i|^2 }{N_0} \right) } {
      \sum_{s_i \in S_1} \exp \left( -\frac{ |r - c s_i|^2 }{N_0} \right) }
      \right)
      \f]
      where \f$s_i \in S_0\f$ denotes a constellation symbol with the i-th bit equal to zero.
    
      \param rx_symbols The received noisy constellation symbols, \f$r\f$ (complex)
      \param channel The channel coefficients (complex), \f$c\f$
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const;

    /*!
      \brief Approximative soft demodulator for a known channel in AWGN
    
      This function is faster and gives almost no performance degradation compared to the 
      demodulate_soft_bits() function. Only include the largest term in nominator and denominator when calculating
      the log-likelihood function.
    */
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const;

    //! Set the size of the signal constellation
    void set_M(int Mary);

    //vec modulate(const svec &symbolnumbers);
    //svec demodulate(const vec &signal);

  protected:
    //! Number of bits per modulation symbol
    int k;
    //! Number of modulation symbols
    int M;
    //! Bit patterns for symbols in order of symbol number
    bmat bitmap; 
    //! Bit pattern in decimal form ordered and the corresponding symbols
    ivec bits2symbols;
    //! A vector containing the modulation symbols
    vec symbols;
    //! The average signal energy of the constallation
    double average_energy;
    //! Scaling factor used to normalize the average energy to 1
    double scaling_factor;
    //! Matrix where row k contains the constellation symbol numbers where bit k is 0/1
    imat S0; 
    //! Matrix where row k contains the constellation symbol numbers where bit k is 0/1
    imat S1; 
  };

  /*! 
    \ingroup modulators
    \brief QPSK-modulator class.

    Symbol numbering is counter clockwise starting with  \f$(1,1)/\sqrt{2}\f$ as symbol 0.
    The bits are Gray coded onto symbols.
    The energy is normalized to one.

    Example of use:
    \code
    QPSK qpsk;
    bvec bits = "0 0 0 1 1 0 1 1";
    cvec symbols = qpsk.modulate_bits(bits);
    \endcode

    This class can also perform soft demodulation, calculating the log-MAP estimate of the individual bits.
    To use the soft demodulate member functions the
    received symbols shall equal
    \f[ r_k = c_k \times s_k + n_k, \f]
    where \f$c_k\f$  is the complex channel
    gain, \f$s_k\f$ is the transmitted QPSK symbols, and \f$n_k\f$ is the AWGN of the channel (with 
    variance \f$N_0/2\f$ in both the real and the imaginary valued components).
  
    The input samples to the soft demodulate functions should be \f$r_k\f$. It is also assumed that the
    channel estimates are perfect when calculating the soft bits.
  
    When these member functions are used together with MAP-based turbo decoding algoritms then the channel
    reliability factor \f$L_c\f$ of the turbo decoder shall be set to 1. The output from these member
    functions can also be used by a Viterbi decoder. 
  */
  class QPSK : virtual public Modulator {
  public:
    //! Class Constructor
    QPSK() {}
    //! Destructor
    virtual ~QPSK() {}
    //! Returns number of bits per symbol (2 for QPSK)
    virtual double bits_per_symbol() const {return 2;}
    //! Modulation of bits
    virtual void modulate_bits(const bvec &bits, cvec &out) const;
    //! Modulation of bits
    virtual cvec modulate_bits(const bvec &bits) const;
    //! Demodulation of bits
    virtual void demodulate_bits(const cvec &signal, bvec &out) const;
    //! Demodulation of bits
    virtual bvec demodulate_bits(const cvec &signal) const;
    /*! 
      \brief Soft demodulator for AWGN channels
    
      This function calculates the log-MAP estimates assuming equally likely bits transmitted
      \f[
      \log \left( \frac{\Pr(b_0=0|r)}{\Pr(b_0=1|r)} \right)
      = \frac{4 \Re \{r\} }{N_0}
      \f]
      \f[
      \log \left( \frac{\Pr(b_1=0|r)}{\Pr(b_1=1|r)} \right)
      = \frac{4 \Im \{r\} }{N_0}
      \f]
      It is assumed that what is received is \f$r = s + n\f$. \f$s\f$ is the QPSK symbol
      and the mapping between symbols and bits is Gray-coded.
      \param rx_symbols The received noisy constellation symbols, \f$r\f$ (real)
      \param N0 The single sided spectral density of the AWGN noise, \f$n\f$
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits) const;

    /*! 
      \brief Soft demodulator for a known channel in AWGN
    
      This function calculates the log-MAP estimate assuming equally likely bits transmitted
      \f[
      \log \left( \frac{\Pr(b_0=0|r)}{\Pr(b_0=1|r)} \right)
      = \frac{4 \Re \{r c^{*} \} }{N_0}
      \f]
      \f[
      \log \left( \frac{\Pr(b_1=0|r)}{\Pr(b_1=1|r)} \right)
      = \frac{4 \Im \{r c^{*} \} }{N_0}
      \f]
      It is assumed that what is received is the complex-valued model: \f$r = c s + n\f$.
      \param rx_symbols The received noisy constellation symbols, \f$r\f$ (complex)
      \param channel The channel coefficients, \f$c\f$ (complex)
      \param N0 The single sided spectral density of the AWGN noise, \f$n\f$
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const;
  
    /*!
      \brief Soft demodulation. Same as exact soft demodulation give above
    */
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits) const;

    /*!
      \brief Soft demodulation. Same as exact soft demodulation give above
    */
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const;

    //cvec modulate(const svec &symbolnumbers);
    //svec demodulate(const cvec &signal);

  };

  /*! 
    \ingroup modulators
    \brief Mary-PSK modulator.

    M-ary PSK modulator with \f$M = 2^k, \_ k = 1, 2, \ldots \f$.
    Symbol numbering counter clockwise starting from the real axis.
    The bit map is Gray encoded. The symbol energy is 1.

    This class can also perform soft demodulation, calculating the log-MAP estimate of the individual bits.
    To use the soft demodulate member functions the
    received symbols shall equal
    \f[ r_k = c_k \times s_k + n_k, \f]
    where \f$c_k\f$  is the complex channel
    gain, \f$s_k\f$ is the transmitted M-PSK symbols, and \f$n_k\f$ is the AWGN of the channel (with 
    variance \f$N_0/2\f$ in both the real and the imaginary valued components).
  
    The input samples to the soft demodulate functions should be \f$r_k\f$. It is also assumed that the
    channel estimates are perfect when calculating the soft bits.
  
    When these member functions are used together with MAP-based turbo decoding algoritms then the channel
    reliability factor \f$L_c\f$ of the turbo decoder shall be set to 1. The output from these member
    functions can also be used by a Viterbi decoder. 
  */
  class PSK : virtual public Modulator {
  public:
    //! Class constructor
    PSK(int Mary) { set_M(Mary); }
    //! Destructor
    virtual ~PSK() { }
    //! Returns number of bits per symbol (1 for BPSK)
    virtual double bits_per_symbol() const {return k;}
    //! Modulation of bits
    virtual void modulate_bits(const bvec &bits, cvec &out) const;
    //! Modulation of bits
    virtual cvec modulate_bits(const bvec &bits) const;
    //! Demodulation of bits
    virtual void demodulate_bits(const cvec &signal, bvec &out) const;
    //! Demodulation of bits
    virtual bvec demodulate_bits(const cvec &signal) const;

    /*! 
      \brief Soft demodulator for AWGN channels
    
      This function calculates
      \f[ 
      \log \left( \frac{\Pr(b_i=0|r)}{\Pr(b_i=1|r)} \right) =
      \log \left( \frac{
      \sum_{s_i \in S_0} \exp \left( -\frac{ |r - s_i|^2 }{N_0} \right) } {
      \sum_{s_i \in S_1} \exp \left( -\frac{ |r - s_i|^2 }{N_0} \right) }
      \right)
      \f]
      where \f$s_i \in S_0\f$ denotes a constellation symbol with the i-th bit equal to zero.
      This function can be used on channels where the channel gain is \f$c = 1\f$.
    
      \param rx_symbols The received noisy constellation symbols, \f$r\f$
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits) const;

    /*!
      \brief Approximative soft demodulator for AWGN channel
    
      This function is faster and gives almost no performance degradation compared to the 
      demodulate_soft_bits() function. Only include the largest term in nominator and denominator when calculating
      the log-likelihood function.
    */
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits) const;

    /*! 
      \brief Soft demodulator for a known channel in AWGN
    
      This function calculates
      \f[ 
      \log \left( \frac{\Pr(b_i=0|r)}{\Pr(b_i=1|r)} \right) =
      \log \left( \frac{
      \sum_{s_i \in S_0} \exp \left( -\frac{ |r - c s_i|^2 }{N_0} \right) } {
      \sum_{s_i \in S_1} \exp \left( -\frac{ |r - c s_i|^2 }{N_0} \right) }
      \right)
      \f]
      where \f$s_i \in S_0\f$ denotes a constellation symbol with the i-th bit equal to zero.
    
      \param rx_symbols The received noisy constellation symbols, \f$r\f$
      \param channel The channel coefficients (complex), \f$c\f$
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const;

    /*!
      \brief Approximative soft demodulator for a known channel in AWGN
    
      This function is faster and gives almost no performance degradation compared to the 
      demodulate_soft_bits() function. Only include the largest term in nominator and denominator when calculating
      the log-likelihood function.
    */
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const;


    //! Change the size of the signal constellation
    void set_M(int Mary);

    //cvec modulate(const svec &symbolnumbers);
    //svec demodulate(const cvec &signal);

  protected:
    //! Number of bits per modulation symbol
    int k;
    //! Number of modulation symbols
    int M;
    //! bit patterns for symbols in order of symbol number
    bmat bitmap; 
    //! bit pattern in decimal form ordered and the corresponding symbols
    ivec bits2symbols; 
    //! A vector containing the modulation symbols
    cvec symbols;
    //! The average signal energy of the constallation
    double average_energy;
    //! Matrix where row k contains the constellation symbol numbers where bit k is 0
    imat S0; 
    //! Matrix where row k contains the constellation symbol numbers where bit k is 1
    imat S1; 
  };

  /*! 
    \ingroup modulators
    \brief Modulator class for square lattice Mary-QAM signals.
  
    The size of the signal constellation is \f$M = 2^k,\_ k = 2, 3, \ldots \f$
    The symbol values in each dimension is: \f$ \{ \sqrt{M}-1, \ldots, 3, 1, -1, -3, \ldots,-(\sqrt{M}-1) \} \f$

    <h3>Symbol numbering:</h3>
    <ul>
    <li>Symbol     0: \f$(\sqrt{M}-1)+j(\sqrt{M}-1)\f$</li>
    <li>Symbol     1: \f$(\sqrt{M}-3)+j(\sqrt{M}-1)\f$</li>
    <li>...</li>
    <li>Symbol (M-1): \f$-(\sqrt{M}-1)-j(\sqrt{M}-1)\f$.</li>
    </ul>
  
    The symbols are normalized so that the average energy is 1. That is, normalized with \f$ \sqrt{2*(M-1)/3}\f$.

    This class can also perform soft demodulation, calculating the log-MAP estimate of the individual bits.
    To use the soft demodulate member functions the
    received symbols shall equal
    \f[ r_k = c_k \times s_k + n_k, \f]
    where \f$c_k\f$  is the complex channel
    gain, \f$s_k\f$ is the transmitted QAM symbols, and \f$n_k\f$ is the AWGN of the channel (with 
    variance \f$N_0/2\f$ in both the real and the imaginary valued components).
  
    The input samples to the soft demodulate functions should be \f$r_k\f$. It is also assumed that the
    channel estimates are perfect when calculating the soft bits.
  
    When these member functions are used together with MAP-based turbo decoding algoritms then the channel
    reliability factor \f$L_c\f$ of the turbo decoder shall be set to 1. The output from these member
    functions can also be used by a Viterbi decoder. 
  */
  class QAM : virtual public Modulator {
  public:
    //! Class Constructor
    QAM(int Mary) { set_M(Mary); }
    //! Destructor
    virtual ~QAM() { }
    //! Returns number of bits per symbol (1 for BPSK)
    virtual double bits_per_symbol() const {return k;}
    //! Modulation of bits
    virtual void modulate_bits(const bvec &bits, cvec &out) const;
    //! Modulation of bits
    virtual cvec modulate_bits(const bvec &bits) const;
    //! Demodulation of bits
    virtual void demodulate_bits(const cvec &signal, bvec &out) const;
    //! Demodulation of bits
    virtual bvec demodulate_bits(const cvec &signal) const;

    /*! 
      \brief Soft demodulator for AWGN channels
    
      This function calculates
      \f[ 
      \log \left( \frac{\Pr(b_i=0|r)}{\Pr(b_i=1|r)} \right) =
      \log \left( \frac{
      \sum_{s_i \in S_0} \exp \left( -\frac{ |r - s_i|^2 }{N_0} \right) } {
      \sum_{s_i \in S_1} \exp \left( -\frac{ |r - s_i|^2 }{N_0} \right) }
      \right)
      \f]
      where \f$s_i \in S_0\f$ denotes a constellation symbol with the i-th bit equal to zero.
      This function can be used on channels where the channel gain is \f$c = 1\f$.
    
      \param rx_symbols The received noisy constellation symbols, \f$r\f$
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const double N0, vec &soft_bits) const;

    /*!
      \brief Approximative soft demodulator for AWGN channel
    
      This function is faster and gives almost no performance degradation compared to the 
      demodulate_soft_bits() function. Only include the largest term in nominator and denominator when calculating
      the log-likelihood function.
    */
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const double N0, vec &soft_bits) const;

    /*! 
      \brief Soft demodulator for a known channel in AWGN
    
      This function calculates
      \f[ 
      \log \left( \frac{\Pr(b_i=0|r)}{\Pr(b_i=1|r)} \right) =
      \log \left( \frac{
      \sum_{s_i \in S_0} \exp \left( -\frac{ |r - c s_i|^2 }{N_0} \right) } {
      \sum_{s_i \in S_1} \exp \left( -\frac{ |r - c s_i|^2 }{N_0} \right) }
      \right)
      \f]
      where \f$s_i \in S_0\f$ denotes a constellation symbol with the i-th bit equal to zero.
    
      \param rx_symbols The received noisy constellation symbols, \f$r\f$
      \param channel The channel coefficients (complex), \f$c\f$
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const;

    /*!
      \brief Approximative soft demodulator for a known channel in AWGN
    
      This function is faster and gives almost no performance degradation compared to the 
      demodulate_soft_bits() function. Only include the largest term in nominator and denominator when calculating
      the log-likelihood function.
    */
    virtual void demodulate_soft_bits_approx(const cvec &rx_symbols, const cvec &channel, const double N0, vec &soft_bits) const;

    //! Change the size of the signal constellation
    void set_M(int Mary);
    //! Return the constellation symbols used
    cvec get_symbols() const { return symbols; }
    //! Return the bit mapping used in decimal form
    ivec get_bitmap() const { return bits2symbols; }

    //cvec modulate(const svec &symbolnumbers);
    //svec demodulate(const cvec &signal);

  protected:
    //! Number of bits per modulation symbol
    int k;
    //! Number of modulation symbols
    int M;
    //! The square-root of M
    int L;
    //! Bit patterns for symbols in order of symbol number
    bmat bitmap; 
    //! Bit pattern in decimal form ordered and the corresponding symbols
    ivec bits2symbols;
    //! A vector containing the modulation symbols
    cvec symbols;
    //! The average signal energy of the constallation
    double average_energy;
    //! Scaling factor used to normalize the average energy to 1
    double scaling_factor;
    //! Matrix where row k contains the constellation symbol numbers where bit k is 0/1
    imat S0; 
    //! Matrix where row k contains the constellation symbol numbers where bit k is 0/1
    imat S1; 
  };

} // namespace itpp

#endif // #ifndef MODULATOR_H

