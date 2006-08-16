/*!
 * \file 
 * \brief One- and two-dimensional modulators - header file
 * \author Tony Ottosson and Adam Piatyszek
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2006  (see AUTHORS file for a list of contributors)
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

#ifndef MODULATOR_H
#define MODULATOR_H

#include <itpp/base/mat.h>


namespace itpp {

  //! \addtogroup modulators

  /*! 
    \ingroup modulators
    \brief General modulator for one dimensional (1D) signal constellations.

    The Modulator_1D class is designed for modeling any kind of 1D (real
    based) signal constellation. Thus, it can be set to represent 1D
    constellations like BPSK or PAM.

    The constellation of the modulator is described with two vectors. One
    contains real values representing the constellation points, whereas the
    other one includes the corresponding bit sets (in the decimal from).

    \note Only hard demodulation methods are implemented.
  */
  class Modulator_1D {
  public:
    //! Default constructor
    Modulator_1D();
    //! Constructor
    Modulator_1D(const vec& symbols, const ivec& bits2symbols);
    //! Destructor
    virtual ~Modulator_1D() { }

    //! Set the symbol constellation and the corresponding bitmap
    virtual void set(const vec& symbols, const ivec& bits2symbols);

    //! Returns the number of bits per symbol (can be non integral)
    virtual double bits_per_symbol() const { return k; }
    //! Get the symbol constellation
    virtual vec get_symbols() const { return symbols; }
    //! Get the bitmap
    virtual ivec get_bits2symbols() const { return bits2symbols; }

    //! Modulate function for symbols
    virtual vec modulate(const ivec& symbolnumbers) const;
    //! Demodulate function for symbols
    virtual ivec demodulate(const vec& signal) const;

    //! Modulate function for bits
    virtual void modulate_bits(const bvec& bits, vec& out) const;
    //! Modulate function for bits
    virtual vec modulate_bits(const bvec& bits) const;
    //! Demodulate function for bits
    virtual void demodulate_bits(const vec& signal, bvec& out) const;
    //! Demodulate function for bits
    virtual bvec demodulate_bits(const vec& signal) const;

  protected:
    //! Setup indicator
    bool setup_done;
    //! Number of bits per modulation symbol
    int k;
    //! Number of modulation symbols
    int M;
    //! Bit to symbol mapping table (size: M x k)
    bmat bitmap;
    //! Bit to symbol mapping in decimal form (size: M)
    ivec bits2symbols;
    //! Corresponding modulation symbols (size: M)
    vec symbols;
    /*! \brief Matrix where row k contains the constellation points with '0'
      in bit position k */
    imat S0;
    /*! \brief Matrix where row k contains the constellation points with '1'
      in bit position k */
    imat S1;

    //! This function calculates the soft bit mapping matrices S0 and S1
    void calculate_softbit_matricies(const ivec& bits2symbols);
  };


  /*! 
    \ingroup modulators
    \brief General modulator for two-dimensional (2D) signal constellations.

    The Modulator_2D class is designed for modeling any kind of 2D (complex
    based) signal constellation. Therefore it is used as a base class for
    such modulations like PSK or QAM. Moreover, it can be set to represent
    1D constellations like BPSK or PAM (only real part of the complex valued
    signal is used).

    The constellation of the modulator is described with two vectors. One
    contains complex values representing the constellation points, whereas
    the other one includes the corresponding bit sets (in the decimal from).

    Beside hard demapping, this class can also perform soft demodulation. To
    use soft demodulate member functions the received symbols should be
    equal to: \f[r_k = c_k s_k + n_k,\f] where \f$c_k\f$ is the complex
    channel gain, \f$s_k\f$ is the transmitted constellation symbol, and
    \f$n_k\f$ is the AWGN of the channel (with variance \f$N_0/2\f$ in both
    the real and imaginary valued components).
  
    The input samples to the soft demodulate functions should be \f$r_k\f$.
    It is also assumed that the channel estimates are perfect when
    calculating the soft bits.
  
    When these member functions are used together with MAP-based turbo
    decoding algorithms then the channel reliability factor \f$L_c\f$ of the
    turbo decoder shall be set to 1. The output from these member functions 
    can also be used by a Viterbi decoder using an AWGN based metric
    calculation function.
  */
  class Modulator_2D {
  public:
    //! Default constructor
    Modulator_2D();
    //! Constructor
    Modulator_2D(const cvec& symbols, const ivec& bits2symbols);
    //! Destructor
    virtual ~Modulator_2D() {}

    //! Set the constellation to use in the modulator
    virtual void set(const cvec& symbols, const ivec& bits2symbols);
  
    //! Returns number of bits per symbol
    virtual double bits_per_symbol() const { return k; }
    //! Get the symbol values used in the modulator.
    virtual cvec get_symbols() const { return symbols; }
    //! Get the bitmap in the decimal form
    virtual ivec get_bits2symbols() const { return bits2symbols; }

    //! Modulation of symbols
    virtual void modulate(const ivec& symbolnumbers, cvec& out) const;
    //! Modulation of symbols
    virtual cvec modulate(const ivec& symbolnumbers) const;

    //! Demodulation of symbols
    virtual void demodulate(const cvec& signal, ivec& out) const;
    //! Demodulation of symbols
    virtual ivec demodulate(const cvec& signal) const;
  
    //! Modulation of bits
    virtual void modulate_bits(const bvec& bits, cvec& out) const;
    //! Modulation of bits
    virtual cvec modulate_bits(const bvec& bits) const;

    //! Hard demodulation of bits
    virtual void demodulate_bits(const cvec& signal, bvec& bits) const;
    //! Hard demodulation of bits
    virtual bvec demodulate_bits(const cvec& signal) const;
  
    //@{
    /*! 
      \brief Soft demodulator for AWGN channels
    
      This function calculates: 
      \f[\log \left( \frac{P(b_i=0|r)}{P(b_i=1|r)} \right) = \log \left(
      \frac{\sum_{s_i \in S_0} \exp \left( -\frac{|r_k - s_i|^2}{N_0}
      \right)} {\sum_{s_i \in S_1} \exp \left( -\frac{|r_k - s_i|^2}{N_0}
      \right)} \right) \f]

      This function can be used on channels where the channel gain
      \f$c_k=1\f$.
    
      \param rx_symbols The received noisy constellation symbols
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above

      \note For soft demodulation it is suggested to use the N-dimensional
      modulator (Modulator_ND) class instead, which is based on QLLR
      arithmetic and therefore faster and more numerically stable.
    */
    virtual void demodulate_soft_bits(const cvec& rx_symbols, double N0, 
				      vec& soft_bits) const;
    virtual vec demodulate_soft_bits(const cvec& rx_symbols, double N0) const;
    //@}

    //@{
    /*! 
      \brief Soft demodulator for fading channels
    
      This function calculates:
      \f[ \log \left( \frac{P(b_i=0|r)}{P(b_i=1|r)} \right) = \log \left(
      \frac{\sum_{s_i \in S_0} \exp \left( -\frac{|r_k - c_k s_i|^2}{N_0}
      \right)} {\sum_{s_i \in S_1} \exp \left( -\frac{|r_k - c_k
      s_i|^2}{N_0} \right)} \right) \f]

      \param rx_symbols The received noisy constellation symbols \f$r_k\f$
      \param channel The complex valued channel values \f$c_k\f$
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above

      \note For soft demodulation it is suggested to use the N-dimensional
      modulator (Modulator_ND) class instead which is based on QLLR
      arithmetic and therefore faster and more numerically stable.
    */
    virtual void demodulate_soft_bits(const cvec& rx_symbols,
				      const cvec& channel, 
				      double N0, vec& soft_bits) const;
    virtual vec demodulate_soft_bits(const cvec& rx_symbols, 
				     const cvec& channel,
				     double N0) const;
    //@}

    //@{
    /*!
      \brief Approximate soft demodulator for AWGN channels. 
    
      This function is faster and gives almost no performance degradation
      compared to the <tt>demodulate_soft_bits(const cvec& symbols, double
      N0, vec& soft_bits)</tt> function. This function finds for each bit
      the closest constellation points that have a zero and one in the
      corresponding position. Let \f$d_0\f$ denote the distance to the
      closest zero point and \f$d_1\f$ denote the distance to the closest
      one point for the corresponding bit respectively. This algorithm then
      computes \f[\frac{d_1^2 - d_0^2}{N_0}\f]
    
      \param rx_symbols The received noisy constellation symbols
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits_approx(const cvec& rx_symbols, double N0,
					     vec& soft_bits) const;
    virtual vec demodulate_soft_bits_approx(const cvec& rx_symbols, 
					    double N0) const;
    //@}
  
    //@{
    /*! 
      \brief Approximate soft demodulator for fading channels. 
    
      This function is faster and gives almost no performance degradation
      compared to the <tt>demodulate_soft_bits()</tt> function. This
      function finds for each bit the closest constellation points that have
      a zero and one in the corresponding position. Let \f$d_0 = |r_k - c_k
      s_0|\f$ denote the distance to the closest zero point and \f$d_1 =
      |r_k - c_k s_1|\f$ denote the distance to the closest one point for
      the corresponding bit respectively. This algorithm then computes
      \f[\frac{d_1^2 - d_0^2}{N_0}\f]

      \param rx_symbols The received noisy constellation symbols \f$r_k\f$
      \param channel The complex valued channel values \f$c_k\f$
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits_approx(const cvec& rx_symbols,
					     const cvec& channel,
					     double N0, vec& soft_bits) const;
    virtual vec demodulate_soft_bits_approx(const cvec& rx_symbols, 
					    const cvec& channel,
					    double N0) const;
    //@}
  

  protected:
    //! Setup indicator
    bool setup_done;
    //! Number of bits per modulation symbol
    int k;
    //! Number of modulation symbols
    int M;
    //! Bit to symbol mapping table (size: M x k)
    bmat bitmap;
    //! Bit to symbol mapping in decimal form (size: M)
    ivec bits2symbols;
    //! Corresponding modulation symbols (size: M)
    cvec symbols;
    /*! \brief Matrix where row k contains the constellation points with '0'
      in bit position k */
    imat S0;
    /*! \brief Matrix where row k contains the constellation points with '1'
      in bit position k */
    imat S1;

    //! This function calculates the soft bit mapping matrices S0 and S1
    void calculate_softbit_matricies(const ivec& bits2symbols);
  };


  // ----------------------------------------------------------------------
  // QAM : Modulator_2D
  // ----------------------------------------------------------------------

  /*!
    \ingroup modulators
    \brief Modulator class for square lattice M-ary QAM signals.
  
    The size of the QAM constellation is \f$M = 2^k\f$, where \f$k = 1, 2,
    \ldots \f$. Symbol values in each dimension are: \f$\{-(\sqrt{M}-1),
    \ldots, -3, -1, 1, 3, \ldots, (\sqrt{M}-1)\}\f$. The bitmap is Gray
    encoded. Symbols are normalized so that the average energy is 1. That
    is, normalized with \f$\sqrt{2(M-1)/3}\f$.

    Beside hard demapping, this class can also perform soft demodulation
    calculating the log-MAP estimate of the individual bits. To use soft
    demodulate member functions the received symbols should be equal to:
    \f[r_k = c_k s_k + n_k,\f] where \f$c_k\f$ is the complex channel gain,
    \f$s_k\f$ is the transmitted QAM symbol, and \f$n_k\f$ is the AWGN of
    the channel (with variance \f$N_0/2\f$ in both real and imaginary valued
    components).
  
    The input samples to the soft demodulate functions should be \f$r_k\f$.
    It is also assumed that the channel estimates are perfect when
    calculating the soft bits.
  
    When these member functions are used together with MAP-based turbo
    decoding algorithms then the channel reliability factor \f$L_c\f$ of the
    turbo decoder should be set to 1. The output from these member functions
    can also be used by a Viterbi decoder.
  */
  class QAM : public Modulator_2D {
  public:
    //! Class Constructor
    QAM(int M) { set_M(M); }
    //! Destructor
    virtual ~QAM() { }
    //! Change the size of the signal constellation
    void set_M(int M);

    //! Modulation of bits
    void modulate_bits(const bvec& bits, cvec& out) const;
    //! Modulation of bits
    cvec modulate_bits(const bvec& bits) const;

    //! Hard demodulation of bits
    void demodulate_bits(const cvec& signal, bvec& bits) const;
    //! Hard demodulation of bits
    bvec demodulate_bits(const cvec& signal) const;

  protected:
    //! The square-root of M
    int L;
    //! Scaling factor of square QAM constellation (sqrt((M-1)*2/3))
    double scaling_factor;
  };


  // ----------------------------------------------------------------------
  // PSK : Modulator_2D
  // ----------------------------------------------------------------------

  /*! 
    \ingroup modulators
    \brief M-ary PSK modulator.

    This class implements an M-ary PSK modulator with \f$M = 2^k\f$, where
    \f$k = 1, 2, \ldots \f$. The symbol numbering is counter clockwise
    starting from the real axis, i.e. symbol \f$(1, 0)\f$. The bitmap is
    Gray encoded. The symbol energy is normalized to 1.

    Beside hard demapping, this class can also perform soft demodulation,
    calculating the log-MAP estimate of the individual bits. To use the soft
    demodulate member functions the received symbols shall equal \f[ r_k =
    c_k \times s_k + n_k, \f] where \f$c_k\f$ is the complex channel gain,
    \f$s_k\f$ is the transmitted M-PSK symbols, and \f$n_k\f$ is the AWGN of
    the channel (with variance \f$N_0/2\f$ in both the real and the
    imaginary valued components).
  
    The input samples to the soft demodulate functions should be \f$r_k\f$.
    It is also assumed that the channel estimates are perfect when
    calculating the soft bits.
  
    When these member functions are used together with MAP-based turbo
    decoding algorithms then the channel reliability factor \f$L_c\f$ of the
    turbo decoder shall be set to 1. The output from these member functions
    can also be used by a Viterbi decoder.
  */
  class PSK : public Modulator_2D {
  public:
    //! Class constructor
    PSK(int M) { set_M(M); }
    //! Destructor
    virtual ~PSK() { }
    //! Change the size of the signal constellation
    void set_M(int M);

    //! Modulation of bits
    void modulate_bits(const bvec& bits, cvec& out) const;
    //! Modulation of bits
    cvec modulate_bits(const bvec& bits) const;

    //! Hard demodulation of bits
    void demodulate_bits(const cvec& signal, bvec& bits) const;
    //! Hard demodulation of bits
    bvec demodulate_bits(const cvec& signal) const;
  };


  // ----------------------------------------------------------------------
  // QPSK : PSK : Modulator_2D
  // ----------------------------------------------------------------------

  /*! 
    \ingroup modulators
    \brief QPSK-modulator class.

    This is a specialized version of the PSK modulator with \f$M = 4\f$
    constellation points. Symbol numbering is counter clockwise starting
    from the real axis. Bits are Gray coded onto symbols. Symbol energy is
    normalized to 1.

    Example of use:
    \code
    QPSK qpsk;
    bvec bits = "0 0 0 1 1 0 1 1";
    cvec symbols = qpsk.modulate_bits(bits);
    \endcode

    Beside hard demapping, this class can also perform soft demodulation
    calculating the log-MAP estimate of the individual bits. To use the soft
    demodulate member functions the received symbols shall equal \f[r_k =
    c_k s_k + n_k,\f] where \f$c_k\f$ is the complex channel gain, \f$s_k\f$
    is the transmitted QPSK symbol, and \f$n_k\f$ is the AWGN of the channel
    (with variance \f$N_0/2\f$ in both the real and the imaginary valued
    components).
  
    The input samples to the soft demodulate functions should be \f$r_k\f$.
    It is also assumed that the channel estimates are perfect when
    calculating the soft bits.
  
    When these member functions are used together with MAP-based turbo
    decoding algorithms then the channel reliability factor \f$L_c\f$ of the
    turbo decoder shall be set to 1. The output from these member functions
    can also be used by a Viterbi decoder.
  */
  class QPSK : public PSK {
  public:
    //! Class Constructor
    QPSK(): PSK(4) {}
    //! Destructor
    virtual ~QPSK() {}
  };


  // ----------------------------------------------------------------------
  // BPSK : PSK : Modulator_2D
  // ----------------------------------------------------------------------

  /*! 
    \ingroup modulators
    \brief BPSK Modulator Class
  
    This is a specialized version of the PSK modulator with \f$M = 2\f$
    constellation points. Symbols used are \f$\{1, -1\}\f$. 

    Bit mapping:
    - \f$0 \rightarrow 1\f$ 
    - \f$1 \rightarrow -1\f$.
  
    Example of use:
    \code
    BPSK bpsk;
    bvec bits = "1 0 0 1 1 0 1 0 1 0 1 1 1 0";
    vec symbols = bpsk.modulate_bits(bits);
    \endcode

    \note Although the constellation points can be represented in the real
    domain only, this class use complex based interface to be compatible
    with other PSK and QAM based modulators.
  */
  class BPSK : public PSK {
  public:
    //! Constructor
    BPSK(): PSK(2) {}
    //! Destructor
    virtual ~BPSK() {}

    //! Modulate bits into BPSK symbols in real domain
    void modulate_bits(const bvec& bits, vec& out) const;
    //! Demodulate noisy BPSK symbols in real domain into bits
    void demodulate_bits(const vec& signal, bvec& out) const;

    //! Modulate bits into BPSK symbols in complex domain
    void modulate_bits(const bvec& bits, cvec& out) const;
    //! Modulate bits into BPSK symbols  in complex domain
    cvec modulate_bits(const bvec& bits) const;
    //! Demodulate noisy BPSK symbols in complex domain into bits
    void demodulate_bits(const cvec& signal, bvec& out) const;
    //! Demodulate noisy BPSK symbols in complex domain into bits
    bvec demodulate_bits(const cvec& signal) const;

    //@{
    /*! 
      \brief Soft demodulator for AWGN channel
    
      This function calculates the log-MAP estimate assuming equally likely
      bits transmitted:
      \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) = \frac{4 \Re\{r\}}
      {N_0}\f]
      It is assumed that the following signal is received: \f$r = b + n\f$.

      \param rx_symbols The received noisy constellation symbols, \f$r\f$
      (complex but symbols in real part) 
      \param N0 The single sided spectral density of the AWGN noise, \f$n\f$
      \param soft_bits The soft bits calculated using the expression above

      \note For soft demodulation it is suggested to use the N-dimensional
      modulator (Modulator_ND) class instead which is based on QLLR
      arithmetic and therefore faster and more numerically stable.
    */
    virtual void demodulate_soft_bits(const cvec& rx_symbols, double N0,
				      vec& soft_bits) const;
    vec demodulate_soft_bits(const cvec& rx_symbols, double N0) const;

    virtual void demodulate_soft_bits_approx(const cvec& rx_symbols,
					     double N0, vec& soft_bits) const;
    vec demodulate_soft_bits_approx(const cvec& rx_symbols, double N0) const;
    //@}

    //@{
    /*! 
      \brief Soft demodulator for a known channel in AWGN      

      This function calculates the log-MAP estimate assuming equally likely
      bits transmitted: 
      \f[\log \left( \frac{P(b=0|r)}{P(b=1|r)} \right) = 
      \frac{4 \Re\{r c^{*}\}}{N_0}\f]
      It is assumed that the following signal is received:
      \f$r = c b + n\f$.
 
      \param rx_symbols The received noisy constellation symbols, \f$r\f$
      (complex but symbols in real part) 
      \param channel The channel coefficients, \f$c\f$ (complex)
      \param N0 The single sided spectral density of the AWGN noise, \f$n\f$
      \param soft_bits The soft bits calculated using the expression above

      \note For soft demodulation it is suggested to use the N-dimensional
      modulator (Modulator_ND) class instead which is based on QLLR
      arithmetic and therefore faster and more numerically stable.
    */
    virtual void demodulate_soft_bits(const cvec& rx_symbols,
				      const cvec& channel, double N0,
				      vec& soft_bits) const;
    vec demodulate_soft_bits(const cvec& rx_symbols, const cvec& channel,
			     double N0) const;

    virtual void demodulate_soft_bits_approx(const cvec& rx_symbols, 
					     const cvec& channel, double N0,
					     vec& soft_bits) const;
    vec demodulate_soft_bits_approx(const cvec& rx_symbols, const cvec& channel,
				    double N0) const;
    //@}
  };


  // ----------------------------------------------------------------------
  // PAM : Modulator_2D
  // ----------------------------------------------------------------------

  /*! 
    \ingroup modulators
    \brief M-ary PAM modulator
  
    This class implements an M-ary PAM modulator with the following signal
    values: \f$\{-(M-1), \ldots, -3, -1, 1, 3, \ldots, (M-1)\}\f$. Symbol
    numbering is from left to right in the increasing order. Gray encoded
    bitmap is used.

    The symbols are normalized so that the average energy is 1. That is,
    normalized with \f$ \sqrt{(M^2-1)/3}\f$.

    \note Although the constellation points can be represented in the real
    domain only, this class use complex based interface to be compatible
    with other PSK and QAM based modulators.
  */
  class PAM : public Modulator_2D {
  public:
    //! Constructor
    PAM(int M) { set_M(M); }
    //! Destructor
    virtual ~PAM() {}
    //! Set the size of the signal constellation
    void set_M(int M);

    //! Modulate bits into PAM symbols in real domain
    void modulate_bits(const bvec& bits, vec& out) const;
    //! Hard demodulation of PAM symbols in real domain to bits
    void demodulate_bits(const vec& signal, bvec& out) const;

    //! Modulate bits into PAM symbols in complex domain
    void modulate_bits(const bvec& bits, cvec& out) const;
    //! Modulate bits into PAM symbols in complex domain
    cvec modulate_bits(const bvec& bits) const;
    //! Hard demodulation of PAM symbols in complex domain to bits
    void demodulate_bits(const cvec& signal, bvec& out) const;
    //! Hard demodulation of PAM symbols in complex domain to bits
    bvec demodulate_bits(const cvec& signal) const;

    //@{
    /*! 
      \brief Soft demodulator for AWGN channels
    
      This function calculates:
      \f[\log \left( \frac{P(b_i=0|r)}{P(b_i=1|r)} \right) = \log \left(
      \frac{\sum_{s_i \in S_0} \exp \left( -\frac{|r - s_i|^2}{N_0} \right)}
      {\sum_{s_i \in S_1} \exp \left( -\frac{|r - s_i|^2}{N_0} \right)}
      \right) \f]
      This function can be used on channels where the channel gain is
      \f$c=1\f$.
    
      \param rx_symbols The received noisy constellation symbols, \f$r\f$
      (complex, but symbols are real) 
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above

      \note For soft demodulation it is suggested to use the N-dimensional
      modulator (Modulator_ND) class instead which is based on QLLR
      arithmetic and therefore faster and more numerically stable.
    */
    virtual void demodulate_soft_bits(const cvec& rx_symbols, double N0,
				      vec& soft_bits) const;
    virtual vec demodulate_soft_bits(const cvec& rx_symbols, double N0) const;
    //@}

    //@{
    /*!
      \brief Approximate soft demodulator for AWGN channel
    
      This function is faster and gives almost no performance degradation
      compared to the <tt>demodulate_soft_bits()</tt> function. This
      function finds for each bit the closest constellation points that have
      a zero and one in the corresponding position. Let \f$d_0\f$ denote the
      distance to the closest zero point and \f$d_1\f$ denote the distance
      to the closest one point for the corresponding bit respectively. This
      algorithm then computes \f[\frac{d_1^2 - d_0^2}{N_0}\f]
    
      \param rx_symbols The received noisy constellation symbols
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits_approx(const cvec& rx_symbols, double N0,
					     vec& soft_bits) const;
    virtual vec demodulate_soft_bits_approx(const cvec& rx_symbols, 
					    double N0) const;
    //@}


    //@{
    /*! 
      \brief Soft demodulator for a known channel in AWGN
    
      This function calculates:
      \f[ \log \left( \frac{P(b_i=0|r)}{P(b_i=1|r)} \right) = \log \left(
      \frac{\sum_{s_i \in S_0} \exp \left( -\frac{|r - c s_i|^2}{N_0}
      \right)}{\sum_{s_i \in S_1} \exp \left( -\frac{|r - c s_i|^2}{N_0}
      \right)} \right) \f]
    
      \param rx_symbols The received noisy constellation symbols, \f$r\f$
      (complex) 
      \param channel The channel coefficients (complex), \f$c\f$
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above

      \note For soft demodulation it is suggested to use the N-dimensional
      modulator (Modulator_ND) class instead which is based on QLLR
      arithmetic and therefore faster and more numerically stable.
    */
    virtual void demodulate_soft_bits(const cvec& rx_symbols,
				      const cvec& channel, double N0,
				      vec& soft_bits) const;
    virtual vec demodulate_soft_bits(const cvec& rx_symbols, 
				     const cvec& channel, double N0) const;
    //@}

    //@{
    /*!
      \brief Approximate soft demodulator for a known channel in AWGN
    
      This function is faster and gives almost no performance degradation
      compared to the <tt>demodulate_soft_bits()</tt> function. This
      function finds for each bit the closest constellation points that have
      a zero and one in the corresponding position. Let \f$d_0 = |r_k - c_k
      s_0|\f$ denote the distance to the closest zero point and \f$d_1 =
      |r_k - c_k s_1|\f$ denote the distance to the closest one point for
      the corresponding bit respectively. This algorithm then computes
      \f[\frac{d_1^2 - d_0^2}{N_0}\f]

      \param rx_symbols The received noisy constellation symbols \f$r_k\f$
      \param channel The complex valued channel values \f$c_k\f$
      \param N0 The single sided spectral density of the AWGN noise
      \param soft_bits The soft bits calculated using the expression above
    */
    virtual void demodulate_soft_bits_approx(const cvec& rx_symbols,
					     const cvec& channel, double N0,
					     vec& soft_bits) const;
    virtual vec demodulate_soft_bits_approx(const cvec& rx_symbols, 
					    const cvec& channel,
					    double N0) const;
    //@}

  protected:
    //! Scaling factor used to normalize the average energy to 1
    double scaling_factor;
  };



} // namespace itpp

#endif // #ifndef MODULATOR_H
