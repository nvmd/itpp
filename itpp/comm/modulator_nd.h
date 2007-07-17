/*!
 * \file
 * \brief Definition of vector ("MIMO") modulator classes
 * \author  Erik G. Larsson
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

#ifndef MODULATOR_ND_H
#define MODULATOR_ND_H

#include <itpp/base/vec.h>
#include <itpp/comm/llr.h>

namespace itpp {

  /*! \addtogroup modulators
   */

  // ----------------------------------------------------------------------
  // Modulator_ND
  // ----------------------------------------------------------------------

  /*! \ingroup modulators
    \brief Base class for an N-dimensional (ND)
    vector ("MIMO") modulator. See \c ND_UPAM for examples.

    Can also be used for scalar modulation/demodulation as an
    alternative to \c Modulator_1D or \c Modulator_2D. Mixed use of \c
    Modulator_1D/\c Modulator_2D and \c Modulator_ND is <b>not
    advised</b>.
  */
  class Modulator_ND {
  public:
    //! Constructor
    Modulator_ND(LLR_calc_unit llrcalc_in=LLR_calc_unit()) { llrcalc=llrcalc_in; };
    //! Destructor
    ~Modulator_ND() {};

    //! Get number of dimensions
    int get_dim() { return nt; }

    //! Get LLR calculation unit
    LLR_calc_unit get_llrcalc() const { return llrcalc; }

    //! Set LLR calculation unit
    void set_llrcalc(LLR_calc_unit llrcalc_in) { llrcalc=llrcalc_in; };

    //! Get number of bits per modulation symbol
    ivec get_k() { return k; }

    //! Get number of modulation symbols per dimension
    ivec get_M() { return M; }

  protected:
    //! Number of dimensions
    int nt;
    //! LLR calculation unit
    LLR_calc_unit llrcalc;
    //! Number of bits per modulation symbol
    ivec k;
    //! Number of modulation symbols along each dimension
    ivec M;
    //! Bit mapping table (one table per dimension)
    Vec<bmat> bitmap;
    //! Bit pattern in decimal form ordered and the corresponding symbols (one pattern per dimension)
    Vec<ivec> bits2symbols;

    //! Convert LLR to log-probabilities
    QLLRvec probabilities(QLLR l); // some abuse of what QLLR stands for...

    //! Convert LLR to log-probabilities, vector version
    Vec<QLLRvec> probabilities(QLLRvec &l); // some abuse of what QLLR stands for...

    /*! \brief Update LLR (for internal use)

      This function updates the numerator and denominator in the expression
    \f[
     \log \left( \frac{
    \sum_{s: b_k=0}  \exp (-x^2) P(s) }{
    \sum_{s: b_k=1}  \exp (-x^2) P(s) }  \right)
    \f]

    \param logP_apriori vector of a priori probabilities per bit
    \param numerator the logarithm of the numerator in the above expression
    \param denominator the logarithm of the denominator in the above expression
    \param s the symbol vector
    \param x ADD DOCUMENTATION FOR THIS ITEM
    */
    void update_LLR(Vec<QLLRvec> &logP_apriori, QLLRvec &numerator, QLLRvec &denominator, ivec &s, QLLR x);

    /*! \brief Update LLR, for scalar channel (for internal use)

      This function updates the numerator and denominator in the expression
    \f[
     \log \left( \frac{
    \sum_{s: b_k=0}  \exp (-x^2) P(s) }{
    \sum_{s: b_k=1}  \exp (-x^2) P(s) }  \right)
    \f]

    \param logP_apriori vector of a priori probabilities per bit
    \param numerator the logarithm of the numerator in the above expression
    \param denominator the logarithm of the denominator in the above expression
    \param s the symbol vector
    \param scaled_norm ADD DOCUMENTATION FOR THIS ITEM
    \param j ADD DOCUMENTATION FOR THIS ITEM
    */
    void update_LLR(Vec<QLLRvec> &logP_apriori, QLLRvec &numerator, QLLRvec &denominator,
		    int s, QLLR scaled_norm, int j);

  };

  // ----------------------------------------------------------------------
  // Modulator_NRD
  // ----------------------------------------------------------------------

  /*!
    \ingroup modulators
    \brief Base class for N-dimensional vector ("MIMO") channel modulator/demodulators with
    real-valued components.

    This class can be used to perform modulation and demodulation for
    a matrix (MIMO) channel of the form \f[ y=Hx+e \f] where H is a
    channel matrix of dimension \f$n_r\times n_t\f$, \f$y\f$ is a
    received vector of length \f$n_r\f$, \f$x\f$ is a transmitted
    vector of length \f$n_t\f$ and \f$e\f$ is noise.

    The class supports soft-input soft-output demodulation.  It can also be used
    for scalar modulation to take advantage of this feature.

    Complex MIMO channels can be handled by using the \c Modulator_NCD
    class. Alternatively, if the signal constellation is separable in
    I/Q then the complex channel can be first transformed to a real
    channel \f[ G=[H_r, -H_i; H_i, H_r] \f]

    See \c ND_UPAM for examples.
  */
  class Modulator_NRD : public Modulator_ND {
  public:
    //! Constructor
    Modulator_NRD() {};
    //! Destructor
    ~Modulator_NRD() {};

    //! Get modulation symbols per dimension
    Vec<vec> get_symbols() { return symbols; }

    /*! \brief Modulation of bits

    Returns a vector of modulated symbols.

    \param bits vector of the bits to be modulated
    */
    vec modulate_bits(const bvec &bits) const;

    /*! \brief Soft MAP demodulation for multidimensional channel, by
      "brute-force" enumeration of all constellation points.

    This function computes the LLR values

    \f[
    LLR(k) = \log \left( \frac{
    \sum_{s: b_k=0}  \exp \left( -\frac{ |y - Hs|^2 }{2\sigma^2} \right) P(s) }{
    \sum_{s: b_k=1}  \exp \left( -\frac{ |y - Hs|^2 }{2\sigma^2} \right) P(s) }  \right)
    \f]

    without approximations. It is assumed that H, y and
    s are real-valued. Complex-valued channels can be handled via the \c Modulator_NCD class.

    \param  H                channel matrix (m*n)
    \param  y                received vector (m*1)
    \param  sigma2           noise variance, per real dimension
    \param  LLR_apriori      vector of a priori LLR values per bit
    \param  LLR_aposteriori  vector of a posteriori LLR values

    The function performs an exhaustive search over all possible points s in the n-dimensional constellation.
    This is only feasible for relatively small constellations.
    The Jacobian logarithm is used to compute
    the sum-exp expression.

    */
    void map_demod(QLLRvec &LLR_apriori,  QLLRvec &LLR_aposteriori,  double sigma2,  mat &H, vec &y);

    /*! \brief Soft MAP demodulation for parallel  channels without crosstalk.

    This function is equivalent to \c map_demod with \f$H=\mbox{diag}(h)\f$. However, it is
    much faster (the complexity is linear in the number of subchannels).
    */
    void map_demod(QLLRvec &LLR_apriori,  QLLRvec &LLR_aposteriori, double sigma2, vec &h, vec &y);


    //! Print some properties of the MIMO modulator in plain text (mainly to aid debugging)
    friend std::ostream &operator<<(std::ostream &os, const Modulator_NRD &mod);

  protected:

    //! Vector of modulation symbols (along each dimension)
    Vec<vec> symbols;

    /*! \brief For internal use only.

    Update the residual norm \f$|y-Hs|\f$ when moving from one constellation point to an adjacent point

    \param norm The norm to be updated
    \param k The position where s changed
    \param sold Old value of s[k]
    \param snew New value of s[k]
    \param ytH The vector y'H
    \param HtH The Grammian matrix H'H
    \param s The s-vector
    */
    void update_norm(double &norm, int k, int sold, int snew, vec &ytH, mat &HtH, ivec &s);
  };

  /*!
    \relatesalso Modulator_NRD
    \brief Print some properties of the MIMO modulator in plain text (mainly to aid debugging)
  */
  std::ostream &operator<<(std::ostream &os, const Modulator_NRD &mod);

  // ----------------------------------------------------------------------
  // Modulator_NCD
  // ----------------------------------------------------------------------

  /*!  \ingroup modulators
    \brief Base class for vector ("MIMO")
    channel modulator/demodulators with complex valued components.

    This class is equivalent to \c Modulator_NRD except for that all
    quantities are complex-valued.  See \c ND_UPAM for examples.
  */
  class Modulator_NCD : public Modulator_ND {
  public:
    //! Constructor
    Modulator_NCD() {};
    //! Destructor
    ~Modulator_NCD() {};

    //! Get modulation symbols per dimension
    Vec<cvec> get_symbols() { return symbols; }

    //! Modulation of bits
    cvec modulate_bits(const bvec &bits) const;

    /*! \brief Soft MAP demodulation for multidimensional channel, by
      "brute-force" enumeration of all constellation points.

    This function computes the LLR values

    \f[
    LLR(k) = \log \left( \frac{
    \sum_{s: b_k=0}  \exp \left( -\frac{ |y - Hs|^2 }{\sigma^2} \right) P(s)}{
    \sum_{s: b_k=1}  \exp \left( -\frac{ |y - Hs|^2 }{\sigma^2} \right) P(s)}  \right)
    \f]

    without approximations. It is assumed that H, y and s are
    complex-valued.

    \param  H                channel matrix (m*n)
    \param  y                received vector (m*1)
    \param  sigma2           noise variance, per complex dimension
    \param  LLR_apriori      vector of a priori LLR values per bit
    \param  LLR_aposteriori  vector of a posteriori LLR values

    The function performs an exhaustive search over all possible
    points s in the n-dimensional constellation.  This is only
    feasible for relatively small constellations.  The Jacobian
    logarithm is used to compute the sum-exp expression.

    */
    void map_demod(QLLRvec &LLR_apriori,  QLLRvec &LLR_aposteriori,  double sigma2,  cmat &H, cvec &y) ;

    /*! \brief Soft MAP demodulation for diagonal channels (without crosstalk).

    This function is equivalent to \c map_demod with a diagonal H, but much faster.

    */
    void map_demod(QLLRvec &LLR_apriori,  QLLRvec &LLR_aposteriori, double sigma2, cvec &H, cvec &y) ;

    //! Print some properties of the MIMO modulator in plain text (mainly to aid debugging)
    friend std::ostream &operator<<(std::ostream &os, const Modulator_NCD &mod);

  protected:

    //! Vector of modulation symbols (along each dimension)
    Vec<cvec> symbols;

    //! ADD DOCUMENTATION FOR THIS ITEM
    void update_norm(double &norm, int k, int sold, int snew, cvec &ytH, cmat &HtH, ivec &s);
  };

  /*!
    \relatesalso Modulator_NCD
    \brief Print some properties of the MIMO modulator in plain text (mainly to aid debugging)
  */
  std::ostream &operator<<(std::ostream &os, const Modulator_NCD &mod);

  // ----------------------------------------------------------------------
  // ND_UPAM
  // ----------------------------------------------------------------------

  /*! \brief Multidimensional channel with uniform PAM along each dimension.

  <b>Example: (4*3 matrix channel with 4-PAM)</b>
  \code
  ND_UPAM chan;            // Multidimensional channel with uniform PAM
  chan.set_Gray_PAM(3,4);  // 3-dimensional matrix channel, 4-PAM (2 bits) per dimension
  cout << chan << endl;
  bvec b=randb(3*2);  // 3*2 bits in total
  vec x=chan.modulate_bits(b);
  QLLRvec llr = zeros_i(3*2);
  QLLRvec llr_ap = zeros_i(3*2);  // apriori equiprobable bits
  mat H = randn(4,3);      // 4*3 matrix channel
  double sigma2=0.01; // noise variance
  vec y= H*x + sqrt(sigma2)*randn(4); // add noise
  chan.map_demod(llr_ap, llr, sigma2, H, y);
  cout << "True bits:" << b << endl;
  cout << "LLRs:" << chan.get_llrcalc().to_double(llr) << endl;
  \endcode

  <b>Example: (scalar channel with 8-PAM)</b>
  \code
  ND_UPAM chan;
  chan.set_Gray_PAM(1,8);  // scalar channel, 8-PAM (3 bits)
  cout << chan << endl;
  bvec b=randb(3);
  vec x=chan.modulate_bits(b);
  QLLRvec llr = zeros_i(3);
  QLLRvec llr_ap = zeros_i(3);
  mat H = "1.0";      // scalar channel
  double sigma2=0.01;
  vec y= H*x + sqrt(sigma2)*randn(1); // add noise
  chan.map_demod(llr_ap, llr, sigma2, H, y);
  cout << "True bits:" << b << endl;
  cout << "LLRs:" << chan.get_llrcalc().to_double(llr) << endl;
  \endcode

  \ingroup modulators
  */
  class ND_UPAM : public Modulator_NRD {
  public:
    //! Constructor
    ND_UPAM(int nt_in=1, int Mary=2);
    //! Destructor
    ~ND_UPAM() {};

    //! Set component modulators to M-PAM with Gray mapping
    void set_Gray_PAM(int nt_in=1, int Mary=2);

    //! Set component modulators to M-PAM with Gray mapping, different M per component
    void set_Gray_PAM(int nt_in=1, ivec Mary="2");

    /*! \brief Sphere decoding

    This function solves the integer-constrained minimization problem
    \f[ \mbox{min} |y-Hs| \f] with respect to \f$s\f$ using sphere the
    decoding algorithm and the Schnorr-Eucner search strategy (see
    source code for further implementation notes).  The function
    starts with an initial search radius and increases it with a
    factor ("stepup") until the search succeeds.

    \param  y              received data vector (\f$n_r\times 1\f$)
    \param  H              channel matrix (\f$n_r\times n_t\f$)
    \param  rmax           maximum possible sphere radius to try
    \param  rmin           sphere radius in the first try
    \param  stepup         factor with which the sphere radius is increased
                           if the search fails
    \param  detected_bits  result of the search (hard decisions only, QLLR
                           for a sure "1" is set to 1000)

    The function returns 0 if the search suceeds, and -1 otherwise.
    */
    int sphere_decoding(vec &y, mat &H, double rmin, double rmax,
			double stepup, QLLRvec &detected_bits);

    //! ZF demodulation (NOT IMPLEMENTED YET)
    void ZF_demod(QLLRvec &detected_bits, double sigma2, mat &H, vec &y);

    //! MMSE demodulation (NOT IMPLEMENTED YET)
    void MMSE_demod(QLLRvec &detected_bits, double sigma2, mat &H, vec &y);


  private:
    // Sphere decoding search with Schnorr Eucner strategy.
    int sphere_search_SE(vec &y_in, mat &H, imat &zrange, double r, ivec &zhat);

    vec spacing;  // spacing between the constellation points

    inline int sign_i(int a) { return (a>0 ? 1 : -1); };
    inline int sign_i(double a) { return (a>0 ? 1 : -1); };
  };

  // ----------------------------------------------------------------------
  // ND_UQAM
  // ----------------------------------------------------------------------

  /*! \brief Complex MIMO with uniform QAM  per dimension

  \ingroup modulators
   */
  class ND_UQAM : public Modulator_NCD {
  public:
    //! Constructor
    ND_UQAM(int nt_in=1, int Mary=4);
    //! Destructor
    ~ND_UQAM() {};

    //! Set component modulators to M-QAM with Gray mapping
    void set_Gray_QAM(int nt_in=1, int Mary=4);

    //! Set component modulators to M-QAM with Gray mapping, different M per component
    void set_Gray_QAM(int nt_in=1, ivec Mary="4");

  protected:
    ivec L;  //!< the square root of M

  };

  // ----------------------------------------------------------------------
  // ND_UPSK
  // ----------------------------------------------------------------------

  /*! Complex MIMO with uniform PSK per dimension

  \ingroup modulators
  */
  class ND_UPSK : public Modulator_NCD {
  public:
    //! Constructor
    ND_UPSK(int nt_in=1, int Mary=2);
    //! Destructor
    ~ND_UPSK() {};

    //! Set component modulators to M-QAM with Gray mapping
    void set_Gray_PSK(int nt_in=1, int Mary=4);

    //! Set component modulators to M-QAM with Gray mapping, different M per component
    void set_Gray_PSK(int nt_in=1, ivec Mary="4");
  };


} // namespace itpp

#endif // #ifndef MODULATOR_ND_H

