/*!
 * \file 
 * \brief Implementation of Communication Channel classes and functions
 * \author Tony Ottosson, Pal Frenger and Zbigniew Dlugaszewski
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

#include <itpp/comm/channel.h>
#include <itpp/base/binary.h>
#include <itpp/base/stat.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/scalfunc.h>
#include <itpp/base/elmatfunc.h>
#include <itpp/base/specmat.h>
#include <itpp/base/filter.h>
#include <itpp/base/transforms.h>
#include <itpp/base/bessel.h>

namespace itpp {

  vec jake_filter(double NormFDopp, int order)
  {
    int i, L = (int)std::floor(double(order)/2);
    double t, x0;
    vec x_pos(L), x_neg(L), x(2*L+1), h(2*L+1);
    for (i=1; i<L+1; i++) {	
      t = double(i);
      x_pos(i-1) = besselj(0.25, 2*pi*NormFDopp*t) / std::sqrt(std::sqrt(t));
    }
    x0 = double(1.468813) * std::sqrt(std::sqrt(NormFDopp));
    x_neg = reverse(x_pos);
    x = concat( concat(x_neg, x0), x_pos);
    h = elem_mult(hamming(2*L+1), x);
    h /= norm(h);
    return h;
  }



  // ------------------------------------------------------------------------------------------------------------------
  Fading_Generator::Fading_Generator(const double norm_doppler, const DOPPLER_SPECTRUM spectrum)
  {
    set_norm_doppler(norm_doppler);
    set_doppler_spectrum(spectrum);
    los_power = 0.0; // no LOS component
    los_dopp = 0.0;
  }

  void Fading_Generator::set_norm_doppler(const double norm_doppler)
  {
    it_assert((norm_doppler >= 0) && (norm_doppler < 1.0), "Fading_Generator: Normalized Doppler must be >= 0 and < 1.");
    n_dopp = norm_doppler; 
    init_flag = false;
  }
  
  void Fading_Generator::set_doppler_spectrum(const DOPPLER_SPECTRUM spectrum)
  {
    dopp_spectrum = spectrum;
    if (spectrum != Rice) { // no LOS component
      los_power = 0.0;
      los_dopp = 0.0;
    }
    init_flag = false;
  }

  void Fading_Generator::set_LOS(const double relative_power, const double relative_doppler)
  {
    it_assert((relative_doppler >= 0) && (relative_doppler <= 1.0), "Relative Doppler must be >=0 and <=1.");
    it_assert(relative_power >= 0.0, "Rice factor need to be >= 0.0");
    it_assert(dopp_spectrum == Rice, "Can only set LOS component if Rice spectrum");
    los_power = relative_power;
    los_dopp = relative_doppler;
    init_flag = false;
  }

  cvec Fading_Generator::generate(const int no_samples)
  {
    cvec output;
    this->generate(no_samples, output);
    return output;
  }

  cvec Fading_Generator::generate(const int no_samples, const int upsampling_factor)
  {
    cvec output;
    this->generate(no_samples, upsampling_factor, output);
    return output;
  }
    
  void Fading_Generator::shift_time_offset(const int no_samples) 
  {
    time_offset += no_samples;
  }

  void Fading_Generator::generate_zero_doppler(const int no_samples, cvec &output)
  {
    output = randn_c(no_samples);
    if(los_power > 0.0) {
      double diffuse = std::sqrt(1.0/(1.0+los_power));
      double direct = diffuse*std::sqrt(los_power);
      for (int i=0; i<no_samples; i++)
	output(i) = diffuse*output(i) + direct*std::complex<double>(std::cos(2*pi*los_dopp*n_dopp*(i+time_offset)),std::sin(2*pi*los_dopp*n_dopp*(i+time_offset)));
    }
    time_offset += no_samples; 
  }

  void Fading_Generator::generate_zero_doppler(const int no_samples, const int upsampling_factor, cvec &output)
  {
    output = randn_c(no_samples);
    if (los_power > 0.0) {
      double diffuse = std::sqrt(1.0/(1.0+los_power));
      double direct = diffuse*std::sqrt(los_power);
      for (int i=0; i<no_samples; i++)
	output(i) = diffuse*output(i) + direct*std::complex<double>(std::cos(2*pi*los_dopp*n_dopp*(i*upsampling_factor+time_offset)),std::sin(2*pi*los_dopp*n_dopp*(i*upsampling_factor+time_offset)));
    }
    time_offset += no_samples * upsampling_factor; 
}
  

  // ------------------------------------------------------------------------------------------------------------------
  Rice_Fading_Generator::Rice_Fading_Generator(const double norm_doppler, const DOPPLER_SPECTRUM spectrum, const int no_freq, const RICE_METHOD method)
    : Fading_Generator(norm_doppler, spectrum)
  {
    set_no_frequencies(no_freq);
    set_method(method);
  }

  void Rice_Fading_Generator::set_no_frequencies(const int no_freq)
  {
    it_assert(no_freq >= 7, "Rice_Fading_Generator: number of doppler frequencies should be at least 7");
    Ni = no_freq;
    init_flag = false;
  }

  void Rice_Fading_Generator::set_method(const RICE_METHOD method)
  {
    // check if this method works for the given spectrum
    rice_method = method;
    init_flag = false;
  }

  int Rice_Fading_Generator::get_no_frequencies()
  {
    return Ni;
  }

  RICE_METHOD Rice_Fading_Generator::get_method()
  {
    return rice_method;
  }

  void Rice_Fading_Generator::init()
  {
    switch(rice_method) {
    case MEDS: // Method of Exact Doppler Spread (MEDS)
      init_MEDS();
      break;

    default:
      it_error("Rice_Fading_Generator: this method is not implemented");
    };

    time_offset = 0; // state in the process
    init_flag = true; // generator ready to use
  }

  void Rice_Fading_Generator::set_time_offset(const int offset)
  {
    it_assert(offset >= 0, "Rice_Fading_Generator: time offset need to be >=0");
    time_offset = offset;
  }
  
  double Rice_Fading_Generator::get_time_offset()
  {
    return time_offset;
  }

  void Rice_Fading_Generator::generate(const int no_samples, cvec &output)
  {
    if (init_flag == false)
      init();

    if (n_dopp == 0.0)
      generate_zero_doppler(no_samples, output);
    else {
      output.set_size(no_samples, false);

      for (int i=0; i<no_samples; i++) {
	output(i) = std::complex<double>( sum( elem_mult( c1, cos(2*pi*f1*n_dopp*(i+time_offset)+th1) ) ), 
				     sum( elem_mult( c2, cos(2*pi*f2*n_dopp*(i+time_offset)+th2) ) ) );
      }
      
      if(los_power > 0.0) { // LOS component exist
	double diffuse = std::sqrt(1.0/(1.0+los_power));
	double direct = diffuse*std::sqrt(los_power);
	for (int i=0; i<no_samples; i++)
	  output(i) = diffuse*output(i) + direct*std::complex<double>(std::cos(2*pi*los_dopp*n_dopp*(i+time_offset)),std::sin(2*pi*los_dopp*n_dopp*(i+time_offset)));
      }
      time_offset += no_samples; 
    }
  }


  void Rice_Fading_Generator::generate(const int no_samples, const int upsampling_factor, cvec &output)
  {
    if (init_flag == false)
      init();
    
    if (n_dopp == 0.0) {
      generate_zero_doppler(no_samples, upsampling_factor, output);
    }
    else {
      output.set_size(no_samples, false);
	
      for (int i=0; i<no_samples; i++) {
	output(i) = std::complex<double>( sum( elem_mult( c1, cos(2*pi*f1*n_dopp*(i*upsampling_factor+time_offset)+th1) ) ), 
				     sum( elem_mult( c2, cos(2*pi*f2*n_dopp*(i*upsampling_factor+time_offset)+th2) ) ) );
      }
        
      if(los_power > 0.0) { // LOS component exist
	double diffuse = std::sqrt(1.0/(1.0+los_power));
	double direct = diffuse*std::sqrt(los_power);
	for (int i=0; i<no_samples; i++)
	  output(i) = diffuse*output(i) + direct*std::complex<double>(std::cos(2*pi*los_dopp*n_dopp*(i*upsampling_factor+time_offset)),std::sin(2*pi*los_dopp*n_dopp*(i*upsampling_factor+time_offset)));
      }
      time_offset += no_samples * upsampling_factor; 
    }
  }


  void Rice_Fading_Generator::init_MEDS()
  {
    vec n;

    switch(dopp_spectrum) {
    case Jakes: // Jakes (Clark) spectrum
    case Rice: // Rice also have a Jakes spectrum component
      n = linspace(1,Ni,Ni);
      f1 = sin(pi/(2*Ni)*(n-0.5));
      c1 = std::sqrt(1.0/Ni)*ones(Ni);
      th1 = randu(Ni)*2*pi;

      n = linspace(1,Ni+1,Ni+1);
      f2 = sin(pi/(2*(Ni+1))*(n-0.5));
      c2 = std::sqrt(1.0/(Ni+1))*ones(Ni+1);
      th2 = randu(Ni+1)*2*pi;
      break;

    default:
      it_error("Rice_Fading_Generator: this spectrum is not implemented for the MEDS Rice fading generator");
    };
  }


  // ------------------------------------------------------------------------------------------------------------------

  FIR_Fading_Generator::FIR_Fading_Generator(const double norm_doppler, const DOPPLER_SPECTRUM spectrum, const int filter_length)
    : Fading_Generator(norm_doppler, spectrum)
  {
    set_filter_length(filter_length);
  }


  void FIR_Fading_Generator::set_filter_length(const int filter_length)
  {
    it_assert(filter_length >= 50, "FIR_Fading_Generator: filter length should be at least 50");
    fir_length = filter_length;
    init_flag = false;
  }

  int FIR_Fading_Generator::get_filter_length()
  {
    return fir_length;
  }

  void FIR_Fading_Generator::init()
  {
    double norm_dopp = n_dopp;
    if (n_dopp > 0.0) {
      switch(dopp_spectrum) {
      case Jakes:
      case Rice: // Rice also has a Jakes spectrum component
	upsample_rate = 1; // upsampling rate      
	while (norm_dopp < 0.1) { // calculate a reasonable upsample rate so that normalized doppler is > 0.1
	  norm_dopp *= 2;
	  upsample_rate *= 2;
	}
	fir_filter.set_coeffs(jake_filter(norm_dopp, fir_length));
	break;

      default:
	it_error("FIR_Fading_Generator: doppler spectrum is not implemented");
      };

      // fill filter with dummy data
      cvec dummy = fir_filter(randn_c(fir_length));
      
      left_overs.set_size(0, false);
    }
    init_flag = true; // generator ready to use
  }

  void FIR_Fading_Generator::generate(const int no_samples, cvec &output)
  {
    if (init_flag == false)
      init();

    if (n_dopp == 0.0)
      generate_zero_doppler(no_samples, output);
    else {
      int no_upsamples = ceil_i(double(no_samples-left_overs.size())/double(upsample_rate)) + 1; // calculate number of samples before upsampling
      
      // should make a smarter interpolation here!!!
      lininterp(fir_filter(randn_c(no_upsamples)), upsample_rate, output);
      output = concat(left_overs, output); // add left-overs from previous filtering
      left_overs = output.right(output.size()-no_samples); // save left-overs for next round of filtering
      output.set_size(no_samples, true);

      if(los_power > 0.0) { // LOS component exist
	double diffuse = std::sqrt(1.0/(1.0+los_power));
	double direct = diffuse*std::sqrt(los_power);
	for (int i=0; i<no_samples; i++)
	  output(i) = diffuse*output(i) + direct*std::complex<double>(std::cos(2*pi*los_dopp*n_dopp*(i+time_offset)),std::sin(2*pi*los_dopp*n_dopp*(i+time_offset)));
      }
      time_offset += no_samples; 
    }
  }


  //! is this really correct???
  void FIR_Fading_Generator::generate(const int no_samples, const int upsampling_factor, cvec &output)
  {
    this->generate(no_samples, output);
  }


  // ------------------------------------------------------------------------------------------------------------------
  IFFT_Fading_Generator::IFFT_Fading_Generator(const double norm_doppler, const DOPPLER_SPECTRUM spectrum)
    : Fading_Generator(norm_doppler, spectrum) { }


  void IFFT_Fading_Generator::init()
  {
    init_flag = true; // generator ready to use
  }


  void IFFT_Fading_Generator::generate(const int no_samples, cvec &output)
  {
    if (init_flag == false)
      init();

    if (n_dopp == 0.0)
      generate_zero_doppler(no_samples, output);
    else {
      switch(dopp_spectrum) {

      case Jakes:
      case Rice: // Rice also has a Jakes spectrum component
	generate_Jakes(no_samples, output);
	break;
      default:
	it_error("IFFT_Fading_Generator: doppler spectrum is not implemented");
      };

      if(los_power > 0.0) { // LOS component exist
	double diffuse = std::sqrt(1.0/(1.0+los_power));
	double direct = diffuse*std::sqrt(los_power);
	for (int i=0; i<no_samples; i++)
	  output(i) = diffuse*output(i) + direct*std::complex<double>(std::cos(2*pi*los_dopp*n_dopp*(i+time_offset)),std::sin(2*pi*los_dopp*n_dopp*(i+time_offset)));
      }
      time_offset += no_samples; 
    }
  }
 
  // Is this really correct???
  void IFFT_Fading_Generator::generate(const int no_samples, const int upsampling_factor, cvec &output)
  {
    this->generate(no_samples, output);
  }


  void IFFT_Fading_Generator::generate_Jakes(const int no_samples, cvec &output)
  {
    int Nfft = pow2i(needed_bits(no_samples));
    double df = 1.0/Nfft;
    int noisesamp = (int)std::ceil(n_dopp/df);
    int no_upsample = 1;
		
    while (noisesamp <= 10) { // if too few samples, increase the FFT size
      Nfft *= 2;
      no_upsample *= 2;
      df = 1.0/double(Nfft);
      noisesamp = (int)std::ceil(n_dopp/df);
      it_assert(no_upsample < 128, "IFFT_Fading_Generator: Too low normalized doppler or too small blocks of data. Results in inefficient algorithm with lots of zero-padding");
    }


    vec Fpos = linspace(0,0.5,Nfft/2+1);
    vec F = concat(Fpos, reverse(-Fpos(1,Nfft/2-1)));
    vec S = zeros(Nfft);
	
    for (int i=0; i<F.size(); i++) {
      if (std::abs(F(i)) < n_dopp)
	S(i) = std::sqrt(1.5 / ( pi*n_dopp*std::sqrt(1-std::pow(F(i)/n_dopp,2))));
      else if (std::abs(F(i)) == n_dopp)
	S(i) = 1000000;
    }

    S /= norm(S,2); S *= Nfft;

    cvec x(Nfft);

    // lots of zeros. Not necessary to do the multiplication for all elements!!!
    // S is converted. Not necessary???
    x = ifft( elem_mult(to_cvec(S), concat(randn_c(noisesamp), zeros_c(Nfft-2*noisesamp), randn_c(noisesamp)) ) ); 
    output = x.mid(0, no_samples);   
  }

  // ------------------------------------------------------------------------------------------------------------------
  Channel_Specification::Channel_Specification(const vec &avg_power_dB, const vec &delay_prof)
  {
    set_channel_profile(avg_power_dB, delay_prof);
  }

  Channel_Specification::Channel_Specification(const CHANNEL_PROFILE profile)
  {
    set_channel_profile(profile);
  }

  void Channel_Specification::set_channel_profile(const vec &avg_power_dB, const vec &delay_prof)
  {
    it_assert(min(delay_prof) == 0, "Minimum relative delay must be 0.");
    it_assert(avg_power_dB.size() == delay_prof.size(), "Power and delay vectors must be of equal length!");
    it_assert(delay_prof(0) == 0, "First tap must be at zero delay");


    N_taps = delay_prof.size();
    // taps should be sorted and unique
    for (int i=1; i<N_taps; i++) {
      it_assert(delay_prof(i) > delay_prof(i-1), "Delays should be sorted and unique");
    }

    a_prof_dB = avg_power_dB;
    d_prof = delay_prof;

    // set doppler spectrum to Jakes per default
    tap_doppler_spectrum.set_size(N_taps, false);
    for (int i=0; i<N_taps; i++)
      tap_doppler_spectrum(i) = Jakes;

    discrete = false;
  }

  void Channel_Specification::set_channel_profile(const CHANNEL_PROFILE profile)
  {
    switch (profile) {
      // -------------- ITU Channel models -----------------
    case ITU_Vehicular_A:
      set_channel_profile( vec("0 -1 -9 -10 -15 -20"), vec("0 310 710 1090 1730 2510") * 1e-9 );
      break;

    case ITU_Vehicular_B:
      set_channel_profile( vec("-2.5 0 -12.8 -10 -25.2 -16"), vec("0 300 8900 12900 17100 20000") * 1e-9 );
      break;
      
    case ITU_Pedestrian_A:
      set_channel_profile( vec("0 -9.7 -19.2 -22.8"), vec("0 110 190 410") * 1e-9 );
      break;
 
    case ITU_Pedestrian_B:
      set_channel_profile( vec("0 -0.9 -4.9 -8 -7.8 -23.9"), vec("0 200 800 1200 2300 3700") * 1e-9 );
      break;

      // -------------- COST259 Channel models -----------------
    case COST259_TUx:
      set_channel_profile( vec("-5.7 -7.6 -10.1 -10.2 -10.2 -11.5 -13.4 -16.3 -16.9 -17.1 -17.4 -19 -19 -19.8 -21.5 -21.6 -22.1 -22.6 -23.5 -24.3"), 
			   vec("217 512 514 517 674 882 1230 1287 1311 1349 1533 1535 1622 1818 1836 1884 1943 2048 2140") * 1e-9 );
      break;

    case COST259_RAx:
      set_channel_profile( vec("-5.2 -6.4 -8.4 -9.3 -10 -13.1 -15.3 -18.5 -20.4 -22.4"), vec("0 42 101 129 149 245 312 410 469 528") * 1e-9 );
      set_doppler_spectrum(0, Rice);
      set_LOS( sqr(0.91/0.41), 0.7); // What should the rice factor be??? Not sure from report!
      break;

    case COST259_HTx:
      set_channel_profile( vec("-3.6 -8.9 -10.2 -11.5 -11.8 -12.7 -13.0 -16.2 -17.3 -17.7 -17.6 -22.7 -24.1 -25.8 -25.8 -26.2 -29 -29.9 -30 -30.7"), 
			   vec("356 441 528 546 609 625 842 916 941 15000 16172 16492 16876 16882 16978 17615 17827 17849 18016") * 1e-9 );
      break;

      // -------------- COST207 Channel models -----------------
    case COST207_RA:
      set_channel_profile( vec("0 -2 -10 -20"), vec("0 200 400 600") * 1e-9 );
      set_doppler_spectrum(0, Rice);
      set_LOS( sqr(0.91/0.41), 0.7);
      break;

    case COST207_RA6:
      set_channel_profile( vec("0 -4 -8 -12 -16 -20"), vec("0 100 200 300 400 500") * 1e-9 );
      set_doppler_spectrum(0, Rice);
      set_LOS( sqr(0.91/0.41), 0.7);
      break;

    case COST207_TU:
      set_channel_profile( vec("-3 0 -2 -6 -8 -10"), vec("0 200 600 1600 2400 5000") * 1e-9 );
      set_doppler_spectrum(2, GaussI);
      set_doppler_spectrum(3, GaussI);
      set_doppler_spectrum(4, GaussII);
      set_doppler_spectrum(5, GaussII);
      break;

    case COST207_TU6alt:
      set_channel_profile( vec("-3 0 -2 -6 -8 -10"), vec("0 200 500 1600 2300 5000") * 1e-9 );
      set_doppler_spectrum(3, GaussI);
      set_doppler_spectrum(4, GaussII);
      set_doppler_spectrum(5, GaussII);
      break;

    case COST207_TU12:
      set_channel_profile( vec("-4 -3 0 -2 -3 -5 -7 -5 -6 -9 -11 -10"), vec("0 200 400 600 800 1200 1400 1800 2400 3000 3200 5000") * 1e-9 );
      set_doppler_spectrum(3, GaussI);
      set_doppler_spectrum(4, GaussI);
      set_doppler_spectrum(5, GaussI);
      set_doppler_spectrum(6, GaussI);
      set_doppler_spectrum(7, GaussI);
      set_doppler_spectrum(8, GaussII);
      set_doppler_spectrum(9, GaussII);
      set_doppler_spectrum(10, GaussII);
      set_doppler_spectrum(11, GaussII);
      break;

    case COST207_TU12alt:
      set_channel_profile( vec("-4 -3 0 -2.6 -3 -5 -7 -5 -6.5 -8.6 -11 -10"), vec("0 200 400 600 800 1200 1400 1800 2400 3000 3200 5000") * 1e-9 );
      set_doppler_spectrum(4, GaussI);
      set_doppler_spectrum(5, GaussI);
      set_doppler_spectrum(6, GaussI);
      set_doppler_spectrum(7, GaussI);
      set_doppler_spectrum(8, GaussII);
      set_doppler_spectrum(9, GaussII);
      set_doppler_spectrum(10, GaussII);
      set_doppler_spectrum(11, GaussII);
      break;

    case COST207_BU:
      set_channel_profile( vec("-3 0 -3 -5 -2 -4"), vec("0 400 1000 1600 5000 6600") * 1e-9 );
      set_doppler_spectrum(2, GaussI);
      set_doppler_spectrum(3, GaussI);
      set_doppler_spectrum(4, GaussII);
      set_doppler_spectrum(5, GaussII);
      break;

    case COST207_BU6alt:
      set_channel_profile( vec("-2.5 0 -3 -5 -2 -4"), vec("0 300 1000 1600 5000 6600") * 1e-9 );
      set_doppler_spectrum(2, GaussI);
      set_doppler_spectrum(3, GaussI);
      set_doppler_spectrum(4, GaussII);
      set_doppler_spectrum(5, GaussII);
      break;

    case COST207_BU12:
      set_channel_profile( vec("-7 -3 -1 0 -2 -6 -7 -1 -2 -7 -10 -15"), vec("0 200 400 800 1600 2200 3200 5000 6000 7200 8200 10000") * 1e-9 );
      set_doppler_spectrum(3, GaussI);
      set_doppler_spectrum(4, GaussI);
      set_doppler_spectrum(5, GaussII);
      set_doppler_spectrum(6, GaussII);
      set_doppler_spectrum(7, GaussII);
      set_doppler_spectrum(8, GaussII);
      set_doppler_spectrum(9, GaussII);
      set_doppler_spectrum(10, GaussII);
      set_doppler_spectrum(11, GaussII);
      break;

    case COST207_BU12alt:
      set_channel_profile( vec("-7.7 -3.4 -1.3 0 -2.3 -5.6 -7.4 -1.4 -1.6 -6.7 -9.8 -15.1"), ivec("0 100 300 700 1600 2200 3100 5000 6000 7200 8100 10000") * 1e-9 );
      set_doppler_spectrum(3, GaussI);
      set_doppler_spectrum(4, GaussI);
      set_doppler_spectrum(5, GaussII);
      set_doppler_spectrum(6, GaussII);
      set_doppler_spectrum(7, GaussII);
      set_doppler_spectrum(8, GaussII);
      set_doppler_spectrum(9, GaussII);
      set_doppler_spectrum(10, GaussII);
      set_doppler_spectrum(11, GaussII);
      break;


    case COST207_HT:
      set_channel_profile( vec("0 -2 -4 -7 -6 -12"), vec("0 200 400 600 15000 17200") * 1e-9 );
      set_doppler_spectrum(4, GaussII);
      set_doppler_spectrum(5, GaussII);
      break;

    case COST207_HT6alt:
      set_channel_profile( vec("0 -1.5 -4.5 -7.5 -8 -17.7"), vec("0 100 300 500 15000 17200") * 1e-9 );
      set_doppler_spectrum(4, GaussII);
      set_doppler_spectrum(5, GaussII);
      break;

    case COST207_HT12:
      set_channel_profile( vec("-10 -8 -6 -4 0 0 -4 -8 -9 -10 -12 -14"), vec("0 200 400 600 800 2000 2400 15000 15200 15800 17200 20000") * 1e-9 );
      set_doppler_spectrum(3, GaussI);
      set_doppler_spectrum(4, GaussI);
      set_doppler_spectrum(5, GaussI);
      set_doppler_spectrum(6, GaussII);
      set_doppler_spectrum(7, GaussII);
      set_doppler_spectrum(8, GaussII);
      set_doppler_spectrum(9, GaussII);
      set_doppler_spectrum(10, GaussII);
      set_doppler_spectrum(11, GaussII);
      break;

    case COST207_HT12alt:
      set_channel_profile( vec("-10 -8 -6 -4 0 0 -4 -8 -9 -10 -12 -14"), vec("0 100 300 500 700 1000 1300 15000 15200 15700 17200 20000") * 1e-9 );
      set_doppler_spectrum(4, GaussI);
      set_doppler_spectrum(5, GaussI);
      set_doppler_spectrum(6, GaussI);
      set_doppler_spectrum(7, GaussII);
      set_doppler_spectrum(8, GaussII);
      set_doppler_spectrum(9, GaussII);
      set_doppler_spectrum(10, GaussII);
      set_doppler_spectrum(11, GaussII);
      break;
    };
  }


  void Channel_Specification::set_doppler_spectrum(DOPPLER_SPECTRUM *tap_spectrum)
  {
    for (int i=0; i<N_taps; i++)
      tap_doppler_spectrum(i) = tap_spectrum[i];
  }

  void Channel_Specification::set_doppler_spectrum(const int tap_number, DOPPLER_SPECTRUM tap_spectrum)
  {
    tap_doppler_spectrum(tap_number) = tap_spectrum;
  }

  void Channel_Specification::set_LOS(const double relative_power, const double norm_doppler)
  {
    it_assert(N_taps >= 1, "Cannot set LOS component if not set channel profile");
    it_assert((norm_doppler >= 0) && (norm_doppler <= 1.0), "Normalized Doppler must be >=0 and <=1.");
    it_assert(relative_power >= 0.0, "Rice factor need to be >= 0.0");
    it_assert(tap_doppler_spectrum(0) == Rice, "Can only set LOS component if Rice spectrum");

    los_power = relative_power;
    los_dopp = norm_doppler;
  }

  void Channel_Specification::get_channel_profile(vec &avg_power_dB, vec &delay_prof)
  {
    avg_power_dB = a_prof_dB;
    delay_prof = d_prof;
  }

  vec Channel_Specification::get_avg_power_dB()
  {
    return a_prof_dB;
  }

  vec Channel_Specification::get_delay_prof()
  {     
    return d_prof;
  }

  DOPPLER_SPECTRUM Channel_Specification::get_doppler_spectrum(const int index)
  {
    it_assert((index >= 0) && (index < N_taps), "Channel_Specification: index of of range");
    return tap_doppler_spectrum(index);
  }

  double Channel_Specification::get_LOS_power()
  {
    return los_power;
  }

  double Channel_Specification::get_LOS_doppler()
  {
    return los_dopp;
  }

  double Channel_Specification::calc_mean_excess_delay()
  {
    vec a_prof = inv_dB(a_prof_dB);
    vec delay_prof = d_prof;

    return ( a_prof*delay_prof / sum(a_prof));
  }

  double Channel_Specification::calc_rms_delay_spread()
  {
    vec a_prof = inv_dB(a_prof_dB);
    vec delay_prof = d_prof;


    double a = ( a_prof*delay_prof / sum(a_prof));
    double b = ( a_prof*sqr(delay_prof) / sum(a_prof) );

    return ( std::sqrt(b-a*a) );
  }

  void Channel_Specification::discretize(const double Ts)
  {
    it_assert(N_taps > 0, "Channel_Specification::discretize: no channel profile specified");
    vec delay_prof(N_taps);
    vec power(N_taps);

    int j = 0, no_taps, j_delay = 0;

    vec a_prof = inv_dB(a_prof_dB); // Convert power profile


    it_assert(d_prof(0) == 0, "Channel_Specification: first tap should be at zero delay");
    delay_prof(0) = d_prof(0);
    power(0) = a_prof(0);

    // Taps within ( (j-0.5)Ts,(j+0.5)Ts] are included in the jth tap
    // 
    for(int i=1; i<N_taps; i++) {
      if( d_prof(i) > (j_delay+0.5)*Ts ) {
	while( d_prof(i) > (j_delay+0.5)*Ts ) { j_delay++; } // first skip empty taps
	// create a new tap at (j+1)Ts
	j++;
	delay_prof(j) = j_delay;
	power(j) = a_prof(i);
      } else {
	power(j) += a_prof(i);
      }
    }

    no_taps = j+1; // number of taps found
    power.set_size(no_taps, true);
    delay_prof.set_size(no_taps, true);
    
    // write over existing channel profile with the discretized version
    a_prof_dB = dB(power);
    d_prof = delay_prof * Ts; // still store in seconds
    N_taps = no_taps;

    discrete = true;    
    discrete_Ts = Ts;

  }

  // ------------------------------------------------------------------------------------------------------------------
  TDL_Channel::TDL_Channel(const double norm_doppler, const vec &avg_power_dB, const ivec &delay_prof)
  {
    set_channel_profile(avg_power_dB, delay_prof);
    set_norm_doppler(norm_doppler);
    method = Rice_MEDS; // default generation method
  }

  TDL_Channel::TDL_Channel(Channel_Specification &channel_spec)
  {
    set_channel_profile(channel_spec);
    set_norm_doppler(0.0); // initially set to 0 doppler.

    // set doppler spectrum
    tap_doppler_spectrum.set_size(N_taps, false);
    for (int i=0; i<N_taps; i++) {
      tap_doppler_spectrum(i) = channel_spec.get_doppler_spectrum(i);
    }

    if (tap_doppler_spectrum(0) == Rice) // set LOS component
      set_LOS(channel_spec.get_LOS_power(), channel_spec.get_LOS_doppler());

    method = Rice_MEDS; // default generation method
  }

  TDL_Channel::~TDL_Channel()
  {
    if (fading_gen.size() > 0) { // delete all old generators
      for (int i = 0; i < fading_gen.size(); i++) { 
	if (fading_gen(i) != NULL)
	  delete fading_gen(i);
	fading_gen(i) = NULL;
      }
    }
  }

  void TDL_Channel::set_channel_profile(const vec &avg_power_dB, const ivec &delay_prof)
  {
    it_assert(min(delay_prof) == 0, "Minimum relative delay must be 0.");
    it_assert(avg_power_dB.size() == delay_prof.size(), "Power and delay vectors must be of equal length!");
    it_assert(delay_prof(0) == 0, "First tap must be at zero delay");

    N_taps = delay_prof.size();
    // taps should be sorted and unique
    for (int i=1; i<N_taps; i++) {
      it_assert(delay_prof(i) > delay_prof(i-1), "Delays should be sorted and unique");
    }

    a_prof = pow(10.0,avg_power_dB/20.0); // Convert power profile to ampiltude profile
    a_prof /= norm(a_prof); // Normalize

    d_prof = delay_prof;

    init_flag = false;
  }

  void TDL_Channel::set_norm_doppler(const double norm_doppler)
  {
    it_assert((norm_doppler >= 0) && (norm_doppler <= 1.0), "TDL_Channel: Normalized Doppler must be >=0 and <=1.");
    n_dopp = norm_doppler;
    init_flag = false;
  }

  void TDL_Channel::set_channel_profile_uniform(const int no_taps)
  {
    it_assert(no_taps >= 1, "Minimum number of taps is 1.");

    vec avg_power_dB = zeros(no_taps);
    ivec delay_prof(no_taps);

    for (int i=0; i<no_taps; i++)
      delay_prof(i) = i;

    set_channel_profile(avg_power_dB, delay_prof);
  }

  // not implemented
  void TDL_Channel::set_channel_profile_exponential(const double delay_spread)
  {
    it_assert(delay_spread > 0.0, "Delay spread must be larger than 0.");


    //a_prof /= norm(a_prof); // Normalize
  }


  void TDL_Channel::set_channel_profile(Channel_Specification &channel_spec)
  {
    vec avg_power_dB;
    vec delay_profile;
    
    it_assert(channel_spec.is_discrete() == true, "TDL_Channel: Channel_Specification has not been discretized");

    channel_spec.get_channel_profile(avg_power_dB, delay_profile);
    set_channel_profile(avg_power_dB, to_ivec(round(delay_profile/channel_spec.get_discrete_Ts())) );
  }


  void TDL_Channel::set_doppler_spectrum(DOPPLER_SPECTRUM *tap_spectrum)
  {
    it_assert(N_taps > 0, "TDL_Channel:: set_doppler_spectrum: channel profile not defined yet");
    tap_doppler_spectrum.set_size(N_taps, false);

    for (int i=0; i<N_taps; i++)
      tap_doppler_spectrum(i) = tap_spectrum[i];

    init_flag = false;
  }

  void TDL_Channel::set_generation_method(const FADING_GENERATION_METHOD generation_method)
  {
    method = generation_method;
    init_flag = false;
  }

  void TDL_Channel::set_LOS(const double relative_power, const double norm_doppler)
  {
    it_assert(N_taps >= 1, "Cannot set LOS component if not set channel profile");
    it_assert((norm_doppler >= 0) && (norm_doppler <= 1.0), "Normalized Doppler must be >=0 and <=1.");
    it_assert(relative_power >= 0.0, "Rice factor need to be >= 0.0");
    it_assert(tap_doppler_spectrum(0) == Rice, "Can only set LOS component if Rice spectrum");

    los_power = relative_power;
    los_dopp = norm_doppler;
    init_flag = false;
  }

  void TDL_Channel::get_channel_profile(vec &avg_power_dB, ivec &delay_prof)
  {
    avg_power_dB = 20 * log10(a_prof);
    delay_prof = d_prof;
  }


  vec TDL_Channel::get_avg_power_dB()
  {
    return ( 20 * log10(a_prof) );
  }

  ivec TDL_Channel::get_delay_prof()
  {
    return d_prof;
  }


  double TDL_Channel::get_norm_doppler()
  {
    return n_dopp;
  }

  double TDL_Channel::get_LOS_power()
  {
    return los_power;
  }

  double TDL_Channel::get_LOS_doppler()
  {
    return los_dopp;
  }

  FADING_GENERATION_METHOD TDL_Channel::get_generation_method()
  {
    return method;
  }

  double TDL_Channel::calc_mean_excess_delay()
  {
    return ( sqr(a_prof)*d_prof / sum_sqr(a_prof));
  }

  double TDL_Channel::calc_rms_delay_spread()
  {
    double a = ( sqr(a_prof)*d_prof / sum_sqr(a_prof));
    double b = ( sqr(a_prof)*sqr(to_vec(d_prof)) / sum_sqr(a_prof) );

    return ( std::sqrt(b-a*a) );
  }

  void TDL_Channel::init()
  {
    it_assert(N_taps > 0, "TDL_Channel:: init: channel profile not defined yet");

    if (fading_gen.size() > 0) { // delete all old generators
      for(int i=0; i<fading_gen.size(); i++)
	{ 
	  if (fading_gen(i) != NULL)
	    delete fading_gen(i);

	  fading_gen(i) = NULL;
	}
    }

    fading_gen.set_size(N_taps, false);

    if (tap_doppler_spectrum.size() == 0) { // doppler-spectrum is not set. Assume Jakes
      tap_doppler_spectrum.set_size(N_taps);
      for(int i=0; i<N_taps; i++)
	tap_doppler_spectrum(i) = Jakes;
    }


    // create all generators and set the parameters
    switch (method) {
    case FIR:
      for(int i=0; i<N_taps; i++) { fading_gen(i) = new FIR_Fading_Generator(n_dopp, tap_doppler_spectrum(i)); }
      break;

    case IFFT:
      for(int i=0; i<N_taps; i++) { fading_gen(i) = new IFFT_Fading_Generator(n_dopp, tap_doppler_spectrum(i)); }
      break;

    case Rice_MEDS:
      // Ni= smallest number of doppler frequencies, increase by 2 for each tap to make taps uncorrelated
      for(int i=0; i<N_taps; i++) { fading_gen(i) = new Rice_Fading_Generator(n_dopp, tap_doppler_spectrum(i), 16+i*2, MEDS); }
      break;
      
    default:
      it_error("TDL_Channel::init(): Generation method is not implemented");
    };

    if (tap_doppler_spectrum(0) == Rice && los_power > 0.0) { // set LOS component
      fading_gen(0)->set_LOS(los_power, los_dopp);
    }

    // initialize all fading generators
    for(int i=0; i<N_taps; i++)
      fading_gen(i)->init();

    init_flag = true;
  }

  void TDL_Channel::generate(const int no_samples, Array<cvec> &channel_coeff)
  {
    if(init_flag == false)
      init();

    channel_coeff.set_size(N_taps, false);

    for (int i=0; i<N_taps; i++)
      channel_coeff(i) = a_prof(i) * fading_gen(i)->generate(no_samples);
  }

  void TDL_Channel::generate(const int no_samples, cmat &channel_coeff)
  {
    if(init_flag == false)
      init();

    channel_coeff.set_size(no_samples, N_taps, false);

    for (int i=0; i<N_taps; i++)
      channel_coeff.set_col(i, a_prof(i) * fading_gen(i)->generate(no_samples));
  }


  // no_samples is the actual number of generated samples (with upsampling factor taken into account)
  // time_offset will be increased by no_samples*upsampling_factor
  void TDL_Channel::generate(const int no_samples, const int upsampling_factor, Array<cvec> &channel_coeff)
  {
    if(init_flag == false)
      init();
    
    it_assert(method == Rice_MEDS, "TDL_Channel::generate: use with Rice_MEDS generation metod only"); 
    
    channel_coeff.set_size(N_taps, false);
    
    for (int i=0; i<N_taps; i++) {
      channel_coeff(i) = a_prof(i) * fading_gen(i)->generate(no_samples, upsampling_factor);
    }
  }

  void TDL_Channel::generate(const int no_samples, const int upsampling_factor, cmat &channel_coeff)
  {
    if(init_flag == false)
      init();

    it_assert(method == Rice_MEDS, "TDL_Channel::generate: use with Rice_MEDS generation metod only"); 

    channel_coeff.set_size(no_samples, N_taps, false);
    
    for (int i=0; i<N_taps; i++) {
      channel_coeff.set_col(i, a_prof(i) * fading_gen(i)->generate(no_samples, upsampling_factor));
    }
  }

  void TDL_Channel::shift_time_offset(const int no_samples)
  {
    for (int i = 0; i < N_taps; i++) {
      fading_gen(i)->shift_time_offset(no_samples);
    }
  }


  void TDL_Channel::filter_known_channel(const cvec &input, cvec &output, const Array<cvec> &channel_coeff)
  {
    int maxdelay = max(d_prof);

    output.set_size(input.size()+maxdelay, false);
    output.zeros();

    for (int i=0; i<N_taps; i++)
      output += concat( zeros_c(d_prof(i)), elem_mult(input, channel_coeff(i)), zeros_c(maxdelay-d_prof(i)) );
  }

  void TDL_Channel::filter_known_channel(const cvec &input, cvec &output, const cmat &channel_coeff)
  {
    int maxdelay = max(d_prof);

    output.set_size(input.size()+maxdelay, false);
    output.zeros();

    for (int i=0; i<N_taps; i++)
      output += concat( zeros_c(d_prof(i)), elem_mult(input, channel_coeff.get_col(i)), zeros_c(maxdelay-d_prof(i)) );
  }

  void TDL_Channel::filter(const cvec &input, cvec &output, Array<cvec> &channel_coeff)
  {
    generate(input.size(), channel_coeff);
    filter_known_channel(input, output, channel_coeff);
  }

  void TDL_Channel::filter(const cvec &input, cvec &output, cmat &channel_coeff)
  {
    generate(input.size(), channel_coeff);
    filter_known_channel(input, output, channel_coeff);
  }

  cvec TDL_Channel::filter(const cvec &input, Array<cvec> &channel_coeff)
  {
    cvec output;
    filter(input, output, channel_coeff);
    return output;
  }

  cvec TDL_Channel::filter(const cvec &input, cmat &channel_coeff)
  {
    cvec output;
    filter(input, output, channel_coeff);
    return output;
  }

  void TDL_Channel::filter(const cvec &input, cvec &output)
  {
    Array<cvec> channel_coeff;
    filter(input, output, channel_coeff);
  }

  cvec TDL_Channel::filter(const cvec &input)
  {
    cvec output;
    filter(input, output);
    return output;
  }


  void TDL_Channel::operator()(const cvec &input, cvec &output, Array<cvec> &channel_coeff)
  {
    filter(input, output, channel_coeff);
  }

  void TDL_Channel::operator()(const cvec &input, cvec &output, cmat &channel_coeff)
  {
    filter(input, output, channel_coeff);
  }


  cvec TDL_Channel::operator()(const cvec &input, Array<cvec> &channel_coeff)
  {
    return filter(input, channel_coeff);
  }

  cvec TDL_Channel::operator()(const cvec &input, cmat &channel_coeff)
  {
    return filter(input, channel_coeff);
  }

  cvec TDL_Channel::operator()(const cvec &input)
  {
    return filter(input);
  }


  void TDL_Channel::calc_impulse_response(const Array<cvec> &channel_coeff, Array<cvec> &impulse_response)
  {
    it_assert(init_flag == true, "calc_impulse_response: TDL_Channel is not initialized");
    it_assert(N_taps == channel_coeff.size(), "calc_impulse_response: number of channel taps do not match");

    int no_samples = channel_coeff(0).size();
    it_assert(no_samples > 0, "calc_impulse_response: channel_coeff must contain samples");

    impulse_response.set_size(no_samples);

    for (int i=0; i<no_samples; i++) {
      impulse_response(i).set_size(d_prof(N_taps-1)+1, false);
      impulse_response(i).zeros();

      for (int j=0; j<N_taps; j++)
	impulse_response(i)(d_prof(j)) = channel_coeff(j)(i);

    }
  }
  
  void TDL_Channel::calc_frequency_response(const Array<cvec> &channel_coeff, Array<cvec> &frequency_response, const int fft_size)
  {
    it_assert(init_flag == true, "calc_frequency_response: TDL_Channel is not initialized");
    it_assert(N_taps == channel_coeff.size(), "calc_frequency_response: number of channel taps do not match");

    int no_samples = channel_coeff(0).size();
    it_assert(no_samples > 0, "calc_frequency_response: channel_coeff must contain samples");

    frequency_response.set_size(no_samples);

    it_assert(fft_size > d_prof(N_taps-1), "calc_frequency_response: fft_size must be larger than the maximum delay in samples");
    cvec impulse_response(fft_size);

    for (int i=0; i<no_samples; i++) {
      impulse_response.zeros();

      for (int j=0; j<N_taps; j++)
	impulse_response(d_prof(j)) = channel_coeff(j)(i);

      fft(impulse_response, frequency_response(i));

    }
  }

  void TDL_Channel::calc_frequency_response(const cmat &channel_coeff, cmat &frequency_response, const int fft_size)
  {
    it_assert(init_flag == true, "calc_frequency_response: TDL_Channel is not initialized");
    it_assert(N_taps == channel_coeff.cols(), "calc_frequency_response: number of channel taps do not match");

    int no_samples = channel_coeff.rows();
    it_assert(no_samples > 0, "calc_frequency_response: channel_coeff must contain samples");

    frequency_response.set_size(fft_size, no_samples, false);

    it_assert(fft_size > d_prof(N_taps-1), "calc_frequency_response: fft_size must be larger than the maximum delay in samples");
    cvec impulse_response(fft_size);
    cvec freq;

    for (int i=0; i<no_samples; i++) {
      impulse_response.zeros();

      for (int j=0; j<N_taps; j++)
	impulse_response(d_prof(j)) = channel_coeff(i, j);

      fft(impulse_response, freq);
      frequency_response.set_col(i, freq);
    }
  }

  //------------------------------------------------------------------------------
  // Binary Symetric Channel, BSC.
  //------------------------------------------------------------------------------


  bvec BSC::operator()(const bvec &input)
  {
    int i, length = input.length();
    bvec output(length);
	
    for (i=0; i<length; i++) {
      if (u() <= p) {
	output(i) = input(i) + bin(1);
      } else {
	output(i) = input(i);
      }
    }
    return output;
  }



  cvec AWGN_Channel::operator()(const cvec &input)
  {
    cvec output = input + sigma*randn_c(input.size());
    return output;
  }

  vec AWGN_Channel::operator()(const vec &input)
  {
    vec output = input + sigma*randn(input.size());
    return output;
  }


} // namespace itpp
