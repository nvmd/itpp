/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2004 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Linear prediction functions, and conversion between common representations of linear predictive parameters.
  \author Thomas Eriksson

  Here we define functions for linear predictive analysis and coding. Several functions
  for computation of lpc parameters are given, together with functions to convert between
  various common representations of the lpc polynomial coefficients. The functionality is
  the same as in the MATLAB functions with the same names.
  The short term for the various parameter types are:
  poly - LPC polynomial coefficients
  ac - autocorrelation coefficients
  rc - reflection coefficients
  lar - log area ratios
  is - inverse sine parameters
  lsf - line spectral frequencies

  $Revision$

  $Date$
*/
 
#ifndef __lpcfunc_h
#define __lpcfunc_h

#include "itpp/base/vec.h"

namespace itpp {

  /*! \addtogroup lpc
   */
  //!@{

  //! Returns a chirped version of the input vector
  vec chirp(const vec &a, double factor);
  //! Spectral distortion between two vectors, in dB.
  double sd(const vec &In1, const vec &In2);
  //! Spectral distortion between two vectors, in dB, up to highest frequency highestfreq.
  double sd(const vec &In1, const vec &In2, double highestfreq);
  //! Computes reflection coefficients from autocorrelation, using the Le-Roux-Guegen algorithm.
  vec lerouxguegenrc(const vec &R, int order);
  //! Levinson    - Levinson-Durbin recursion.
  vec levinson(const vec &R2, int order);
  //! Computes the autocorrelation function
  vec autocorr(const vec &x, int order);
  //!    lpc         - Linear Predictive Coefficients using autocorrelation method.
  vec lpc(const vec &x, int order);
  //!    schurrc     - Schur algorithm.
  vec schurrc(const vec &R, int order);
  //!    ac2rc       - Autocorrelation sequence to reflection coefficients conversion. 
  vec ac2rc(const vec &ac);
  //!    ac2poly     - Autocorrelation sequence to prediction polynomial conversion.
  vec ac2poly(const vec &ac);
  //!    is2rc       - Inverse sine parameters to reflection coefficients conversion.
  vec is2rc(const vec &is);
  //!    lar2rc      - Log area ratios to reflection coefficients conversion.
  vec lar2rc(const vec &lar);
  //!    lsf2poly    - Line spectral frequencies to prediction polynomial conversion.
  vec lsf2poly(const vec &lsf);
  //!    poly2ac     - Prediction polynomial to autocorrelation sequence conversion. 
  vec poly2ac(const vec &poly);
  //!    poly2lsf    - Prediction polynomial to line spectral frequencies conversion.
  vec poly2lsf(const vec &poly);
  //!    poly2rc     - Prediction polynomial to reflection coefficients conversion.
  vec poly2rc(const vec &poly);
  //!    poly2cepstrum     - Prediction polynomial to cepstrum conversion.
  vec poly2cepstrum(const vec &a);
  //!    poly2cepstrum     - Prediction polynomial to cepstrum conversion, to the specified order.
  vec poly2cepstrum(const vec &a, int num);
  //!    cepstrum2poly     - Cepstrum to prediction polynomial conversion.
  vec cepstrum2poly(const vec &c);
  //!    rc2ac       - Reflection coefficients to autocorrelation sequence conversion.
  vec rc2ac(const vec &rc);
  //!    rc2is       - Reflection coefficients to inverse sine parameters conversion.
  vec rc2is(const vec &rc);
  //!    rc2lar      - Reflection coefficients to log area ratios conversion.
  vec rc2lar(const vec &rc);
  //!    rc2poly     - Reflection coefficients to prediction polynomial conversion.
  vec rc2poly(const vec &rc);

  //!@}

} //namespace itpp

#endif // __lpcfunc_h
