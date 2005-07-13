/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2002 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Definition of an Orthogonal Frequency Division Multiplex (OFDM) class
  \author Pål Frenger, Anders Persson and Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __ofdm_h
#define __ofdm_h

#include "itpp/base/vec.h"

namespace itpp {

  /*! 
    \ingroup modulators
    \brief Class for modulating and demodulation of OFDM signals using the FFT

    The modulated signal is normalized taking into account the cyclic prefix
  */
  class OFDM {
  public:
    //! Empty constructor
    OFDM(void) { setup_done = false; }
    //! Constructor \a Nfft is the size of the FFT. \a Ncp is the length of the cyclic prefix. \a Nupsample is the upsampling factor (default=1)
    OFDM(int inNfft, int inNcp, int inNupsample=1);
    //! Return the number of carriers
    int no_carriers() {return Nfft;}
    //! Set parameters
    void set_parameters(const int Nfft, const int Ncp, const int inNupsample=1);
    //! Modulate complex data symbols. Length of \c input must be \c Nfft
    cvec modulate(const cvec &input);
    //! Modulate complex data symbols. Length of \c input must be \c Nfft
    void modulate(const cvec &input, cvec &output);
    //! Demodulate to complex valued symbols. Length of \c input must be \c Nfft+Ncp
    cvec demodulate(const cvec &input);
    //! Demodulate to complex valued symbols. Length of \c input must be \c Nfft+Ncp
    void demodulate(const cvec &input, cvec &output);  
  private:
    double norm_factor;
    bool setup_done;
    int Nfft, Ncp, Nupsample;
  };

} //namespace itpp

#endif // __ofdm_h
