/*---------------------------------------------------------------------------*
 *                                   IT++                                    *
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
  \brief Implementation of an Orthogonal Frequency Division Multiplex (OFDM) class.
  \author Pål Frenger, Anders Persson and Tony Ottosson

  $Revision$

  $Date$
*/

#include <itpp/comm/ofdm.h>
#include <itpp/base/specmat.h>
#include <itpp/base/stat.h>
#include <itpp/base/transforms.h>
#include <itpp/base/operators.h>


namespace itpp { 

  OFDM::OFDM(int inNfft, int inNcp, int inNupsample)
  {
    set_parameters(inNfft, inNcp, inNupsample);
  }

  void OFDM::set_parameters(const int inNfft, const int inNcp, const int inNupsample){
    it_assert(inNfft >= 2, "OFDM: Nfft must be >=2.");
    it_assert(inNcp >= 0 && inNcp <= inNfft, "OFDM: Ncp must be >=0 and <=Nfft.");
    it_assert(inNupsample >= 1 && inNupsample <= 100, "OFDM: Ncp must be >=1 and <=100.");
    Nfft = inNfft;
    Ncp = inNcp;  
    Nupsample = inNupsample;
    norm_factor = sqrt(static_cast<double>(Nupsample*Nfft*Nfft) / (Nfft+Ncp));
    setup_done = true;
  }

  void OFDM::modulate(const cvec &input, cvec &output){
    it_assert(setup_done == true, "OFDM::modulate: You must set the length of the FFT and the cyclic prefix!");
    const int N = input.length() / Nfft;
    it_assert(N*Nfft == input.length(), "OFDM::modulate: Length of input vector is not a multiple of Nfft.");

    output.set_length(Nupsample*N*(Nfft+Ncp));
    cvec outtemp(Nfft);

    for (int i = 0; i < N; i++) {
			outtemp = ifft(concat(input.mid(i*Nfft, Nfft/2), zeros_c(Nfft*(Nupsample-1)), 
														input.mid(i*Nfft+Nfft/2, Nfft/2))) * norm_factor;
      output.replace_mid(Nupsample*(Nfft+Ncp)*i, concat(outtemp.right(Nupsample*Ncp), outtemp));
    }
  }

  cvec OFDM::modulate(const cvec &input)
  {
    cvec output;
    modulate(input, output);
    return output;
  }

  void OFDM::demodulate(const cvec& input, cvec &output){
    it_assert(setup_done == true, "OFDM::demodulate: You must set the length of the FFT and the cyclic prefix!");
    const int N = input.length() / (Nfft+Ncp) / Nupsample;
    it_assert(Nupsample*N*(Nfft+Ncp) == input.length(), "OFDM: Length of input vector is not a multiple of Nfft+Ncp.");

    output.set_length(N*Nfft);
    // normalize also taking the energy loss into the cyclic prefix into account
    for (int i = 0; i < N; i++) {
			cvec x = fft(input.mid(Nupsample*(i*(Nfft+Ncp)+Ncp), Nupsample*Nfft));
			output.replace_mid(Nfft*i, concat(x.left(Nfft/2), x.right(Nfft/2)) / norm_factor);
    }
  }

  cvec OFDM::demodulate(const cvec &input){
    cvec output;
    demodulate(input, output);
    return output;
  }

} //namespace itpp
