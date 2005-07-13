/*---------------------------------------------------------------------------*
 *                                   IT++			             *
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
  \brief Definition of Communication Channel classes and functions
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __channel_h
#define __channel_h

#include "itpp/base/vec.h"
#include "itpp/base/mat.h"
#include "itpp/base/random.h"
#include "itpp/base/array.h"
#include "itpp/base/filter.h"

namespace itpp {

  /*!
    \addtogroup channels
    \brief Communication Channel Models
    \author Tony Ottosson

    \section channel_contents Modelling and simulation of communication channels
    <ul>
    <li> \ref channel_intro
    <ul>
    <li> \ref channel_doppler
    <li> \ref channel_freq_sel
    <li> \ref channel_los
    <li> \ref channel_ref
    </ul>
    </ul>

    \section channel_intro Introduction

    When simulating a communication link a model of the channel behaviour is usually needed. These models
    typically consist of three parts: the propagation attenuation, the shadowing (log-normal fading) and the
    multipath fading (small scale fading). In the following we will focus on the small scale (or multipath) fading.

    Multipath fading is the process where the received signal is a sum of many reflections each with different 
    propagation time, phase and attenuation. The sum signal will vary in time if the receiver (or transmitter) moves or
    if some of the reflectors move. We usually refer to this process as a fading process and try to model it as a
    stochastic process. The most common model is the Rayleigh fading model where the process is modeled as a sum of infinitely
    many (in practise it is enough with only a few) received reflections from all angles (uniformly) around the receiver.
    Mathematically we write the receieved signal, \f$r(t)\f$ as
    \f[
    r(t) = a(t) * s(t) ,
    \f]
    where \f$s(t)\f$ is the transmitted signal and \f$a(t)\f$ is the complex channel coefficient (or fading process). If
    this process is modeled as a Rayleigh fading process then \f$a(t)\f$ is a complex Gaussian process and the envelope
    \f$\|a(t)\|\f$ is Rayleigh distributed. 


    \subsection channel_doppler Doppler

    The speed by which the channel changes is decided by the speed of the mobile (transmitter or receiver or both). This movement will
    cause the channel coefficient, \f$a(t)\f$ to be correlated in time (or equivalently in frequency). Different models exist of this
    correlation but the most common is the classical Jakes model where the correlation function is given as
    \f[
    R(\tau) = E[a^*(t) a(t+\tau)] = J_0(2 \pi f_\mathrm{max} \tau) ,
    \f]
    where \f$f_\mathrm{max}\f$ is the maximum doppler frequency given by
    \f[
    f_\mathrm{max} = \frac{v}{\lambda} = \frac{v}{c_0} f_c .
    \f]
    Here \f$c_0\f$ is the speed of light and \f$f_c\f$ is the carrier frequency. Often the maximum doppler frequency is given
    as the normalized doppler \f$f_\mathrm{max} T_s\f$, where \f$T_s\f$ is the sample duration (often the symbol time) of the simulated system.
    Instead of specifing the correlation fuction \f$R(\tau)\f$ we can specify the doppler spectrum (the fourier transform of \f$R(\tau)\f$).


    \subsection channel_freq_sel Frequency-selective channels

    Since \f$a(t)\f$ affects the transmitted signal as a constant scaling factor at a given time
    this channel model is often refered to as flat-fading (or frequency non-selective fading). On the other hand, if time arrivals
    of the reflections are very different (compared to the sample times) we cannot model the received signal only as a scaled
    version of the transmitted signal. Instead we model the channel as frequency-selective but time-invariant (or at least wide-sense stationary) 
    with the impulse response
    \f[
    h(t) = \sum_{k=0}^{N_\mathrm{taps}-1} a_k \exp (-j \theta_k ) \delta(t-\tau_k) ,
    \f]
    where \f$N_\mathrm{taps}\f$ is the number of channel taps, \f$a_k\f$ is the average amplitude at delay \f$\tau_k\f$, and \f$\theta_k\f$ is the channel
    phase of the \f$k\f$th channel tap. The average power profile,  and the delay profiles are defined as:
    \f[
    \mathbf{a} = [a_0, a_1, \ldots, a_{N_\mathrm{taps}-1}]
    \f]
    and
    \f[
    \mathbf{\tau} = [\tau_0, \tau_1, \ldots, \tau_{N_\mathrm{taps}-1}],
    \f]
    respectively. We assume without loss of generality that \f$\tau_0 = 0\f$ and
    \f$\tau_0 < \tau_1 < \ldots < \tau_{N_\mathrm{taps}-1}\f$.
    Now the received signal is simpy a linear filtering (or convolution of the transmitted signal) where, \f$h(t)\f$ is impulse response of the filter.

    In practise, when simulating a communication link, the impulse response \f$h(t)\f$ is sampled with a sample period \f$T_s\f$ that is related
    to the symbol rate of the system of investigation (often 2-8 times higher). Hence, the impulse respone of the channel can now be
    modeled as a time-discrete filter or tapped-delay line (TDL) where the delays are given as \f$\tau_k = d_k T_s\f$, and \f$d_k\f$ are positive
    integers.


    \subsection channel_los Line of Sight (LOS) or Rice fading
    If there is line of sight (LOS) between transmitter and receiver the first component received (shortest delay) will have a static
    component that only depend on the doppler frequency. In practise the difference in time between the first LOS component and
    the first refelcted components is small and hence in a discretized system the first tap is usually modeled as a LOS component
    and a fading component. Such a process is usually called a Rice fading process.

    The LOS component can be expressed as:
    \f[
    \rho \exp(2 \pi f_\rho t + \theta_\rho) ,
    \f]
    where \f$\rho\f$, \f$f_\rho\f$, and \f$\theta_\rho\f$ are the amplitude, doppler frequency and phase of the LOS component, respectively.
    Instead of stating the amplitude itself the ratio of the LOS power and the fading process power (or relative power), often called
    the Rice factor, is often stated. The doppler frequency is limited by the maximum doppler frequency \f$f_\mathrm{max}\f$ and 
    hence typically the doppler of the LOS is expressed relative to its maximum (common is \f$f_\rho = 0.7 f_\mathrm{max}\f$). 
    The phase is usually assumed to be random and can without loss of generality be set to 0 (does not affect the statistics of the process).


    \subsection channel_ref References

    [Pätzold] Matthias Pätzold, Mobile fading channels, Wiley, 2002.

    [Stüber] Gordon L. Stüber, Principles of mobile communication, 2nd. ed., Kluwer, 2001.

    [Rappaport] Theodore S. Rappaport, Wireless communications: principles and practise, Prentice Hall, 1996.
  */



  /*!
    \ingroup channels
    \brief Jakes spectrum filter

    Function that generates the taps in the Jake-filter.
    order is the number of taps in the filter.
    \a NormFdopp is the normalized doppler frequency, i.e. \a NormFDopp = \a Fd * \a Ts, where
    \a Fd is the acctual Doppler frequency and \a Ts is the sampling interval.
    Returns a vector containing the filter taps of the Jake-filter.
  */
  vec jake_filter(double NormFDopp, int order = 100);

  /*! 
    \brief Predefined doppler spectra
    \ingroup channels
  */
  enum DOPPLER_SPECTRUM {Jakes=0,J=0,Classic=0,C=0, GaussI=1,GI=1,G1=1, GaussII=2,GII=2,G2=2, Rice=3,R=3};

  /*!
    \brief Methods calculation of parameters using the Rice fading generation method
    \ingroup channels
  */
  enum RICE_METHOD {MEDS};
  
  /*!
    \brief Fading generation methods.
    \ingroup channels
  */

  enum FADING_GENERATION_METHOD {IFFT, FIR, Rice_MEDS};
  
  /*!
    \brief Predefined channel profiles. Includes settings for doppler spectrum
    \ingroup channels
  */
  enum CHANNEL_PROFILE {ITU_Vehicular_A, ITU_Vehicular_B, ITU_Pedestrian_A, ITU_Pedestrian_B,
			COST207_RA, COST207_RA6, 
			COST207_TU, COST207_TU6alt, COST207_TU12, COST207_TU12alt, 
			COST207_BU, COST207_BU6alt, COST207_BU12, COST207_BU12alt, 
			COST207_HT, COST207_HT6alt, COST207_HT12, COST207_HT12alt,
			COST259_TUx, COST259_RAx, COST259_HTx};

  /*!
    \ingroup channels
    \brief Fading generator class
    \author Tony Ottosson

    Abstract class defining the interface of a single tap fading generator.
    Parameters that define the generator is the normalized doppler and the
    doppler spectrum. Possible values of doppler spectra are:
    <ul>
    <li> \a Jakes, J, Classic, or C: for the classical Jakes spectrum shape</li>
    <li> \a Rice or R: the classical Jakes spectrum and a direct tap. Observe that the LOS parameters also need to be set.</li>
    <li> \a GaussI, GI, Gauss1, or G1</li>
    <li> \a GaussII, GII, Gauss2, or G2</li>
    </ul>
    
    Two types of genators exist in the litterature: the filter method and the Rice method. The filter method filter 
    a complex Gaussian process to achieve a given doppler spectrum while the Rice method approximate the process as
    a sum of sinusoids.

    References:

    [Pätzold] Matthias Pätzold, Mobile fading channels, Wiley, 2002.
  */
  class Fading_Generator {
  public:
    Fading_Generator(const double norm_doppler=0.0, const DOPPLER_SPECTRUM spectrum=Jakes);
    virtual ~Fading_Generator(){ }
 
    //! Set normalized doppler rate
    void set_norm_doppler(const double norm_doppler);
    //! Set doppler spectrum
    void set_doppler_spectrum(const DOPPLER_SPECTRUM spectrum);
    //! Set LOS component. Used e.g. for Rice doppler spectrum. Rice factor and relative doppler (related to maximum doppler).
    void set_LOS(const double relative_power, const double relative_doppler);

    //! Return normalized doppler rate
    double get_norm_doppler() { return n_dopp; }
    //! Return doppler spectrum
    DOPPLER_SPECTRUM get_doppler_spectrum() { return dopp_spectrum; }
    //! Get relative power of LOS component (Rice factor)
    double get_LOS_power() { return los_power; }
    //! Get relative doppler (compared to maximum doppler) for the LOS component
    double get_LOS_doppler() { return los_dopp; }

    //! Initialize the generator (has to be done before calling generate())
    virtual void init() = 0;

    //! Generate \c no_samples values from the fading process
    virtual void generate(const int no_samples, cvec &output) = 0;
    //! Generate \c no_samples values from the fading process
    cvec generate(const int no_samples);

    //! Generate \c no_samples values from the fading process upsampled by \c upsampling_factor
    virtual void generate(const int no_samples, const int upsampling_factor, cvec &output) = 0;
    //! Generate \c no_samples values from the fading process upsampled by \c upsampling_factor
    cvec generate(const int no_samples, const int upsampling_factor);
   
    //! Shift generator time offset by a number of samples
    void shift_time_offset(const int no_samples);

  protected:
    DOPPLER_SPECTRUM dopp_spectrum; // doppler spectrum
    double n_dopp; // normalized doppler
    double los_dopp; // relative doppler on LOS-component
    double los_power; // Relative power of LOS component compared to diffuse component = Rice factor < 1
    bool init_flag; // signals if generator is initialized or not
    double time_offset; // time offset in samples (the time state in the generator, not used in the filter method other than for LOS)
    void generate_zero_doppler(const int no_samples, cvec &output); // generate fading for n_dopp == 0.0
    void generate_zero_doppler(const int no_samples, const int upsampling_factor, cvec &output); // generate fading for n_dopp == 0.0
  };

  /*!
    \ingroup channels
    \brief Rice type Fading generator class
    \author Tony Ottosson

    A Rice generator is a generator of the form
    \f[
    \tilde \mu_i(t) = \sum_{n=1}^{N_i} c_{i,n} \cos(2\pi f_{i,n} t + \theta_{i,n})
    \f]
    Here \f$c_{i,n}\f$, \f$f_{i,n}\f$, and \f$\theta_{i,n}\f$ are the doppler coefficients, discrete doppler frequencies, and
    doppler phases, respectively. Rice showed that a generator of this form can perfectly model a Gaussian process
    when \f$N_i \leftarrow \infty\f$. When generating a fading pattern we need a complex-valued generator
    \f[
    \tilde \mu(t) = \tilde \mu_1(t) + j \tilde \mu_2(t)
    \f]
    Parameters that define the generator is the normalized doppler and the
    doppler spectrum. Possible values of doppler spectra are:
    <ul>
    <li> \a Jakes, J, Classic, or C: for the classical Jakes spectrum shape</li>
    <li> \a Rice or R: the classical Jakes spectrum and a direct tap. Observe that the LOS parameters also need to be set.</li>
    <li> \a GaussI, GI, Gauss1, or G1</li>
    <li> \a GaussII, GII, Gauss2, or G2</li>
    </ul>
    
    Furthermore also the number of doppler frequencies, \f$N_i\f$ and the method used to calculate the parameters
    \f$c_{i,n}\f$, \f$f_{i,n}\f$, and \f$\theta_{i,n}\f$ are parameters. For now the only method defined for calculating
    the parameters is the Method of Exact Doppler Spread (MEDS). See [Pätzold] for more details.

    References:

    [Pätzold] Matthias Pätzold, Mobile fading channels, Wiley, 2002.
  */
  class Rice_Fading_Generator : public Fading_Generator {
  public:
    //! Set normalized dopper, doppler spectrum, number of doppler frequencies and calculation method
    Rice_Fading_Generator(const double norm_doppler = 0.0, const DOPPLER_SPECTRUM spectrum = Jakes, 
			  const int no_freq = 16, const RICE_METHOD method=MEDS);
    //! Destructor
    virtual ~Rice_Fading_Generator(){ }

    //! Set number of doppler frequencies
    void set_no_frequencies(const int no_freq);
    //! Set calculation method for calculation of doppler frequencies and amplitudes
    void set_method(const RICE_METHOD method);
    //! Get number of doppler frequencies
    int get_no_frequencies();
    //! Get calculation method for calculation of doppler frequencies and amplitudes
    RICE_METHOD get_method();

    //! Initialize the generator (is not needed) and set time offset to 0
    virtual void init();

    //! set time offset in samples
    void set_time_offset(const int offset);
    //! get time offset in samples
    double get_time_offset();

    //! Generate \c no_samples values from the fading process
    virtual void generate(const int no_samples, cvec &output);
    //! Generate \c no_samples values from the fading process upsampled by \c upsampling_factor
    virtual void generate(const int no_samples, const int upsampling_factor, cvec &output);

    //! Generate \c no_samples values from the fading process
    //virtual cvec generate(const int no_samples);

  protected:
    int Ni;
    RICE_METHOD rice_method;
    vec f1, f2, c1, c2, th1, th2; // doppler frequencies, amplitudes and phases (need to be initialized)

    //! Init function for MEDS method
    void init_MEDS();
  };

  /*!
    \ingroup channels
    \brief FIR type Fading generator class
    \author Tony Ottosson

    A FIR generator is a linear finite impulse response (FIR) filter implementation of filter method for generation of a fading process.
    Parameters that define the generator is the normalized doppler and the
    doppler spectrum. Possible values of doppler spectra are:
    <ul>
    <li> \a Jakes, J, Classic, or C: for the classical Jakes spectrum shape</li>
    <li> \a Rice or R: the classical Jakes spectrum and a direct tap. Observe that the LOS parameters also need to be set.</li>
    <li> \a GaussI, GI, Gauss1, or G1</li>
    <li> \a GaussII, GII, Gauss2, or G2</li>
    </ul>
    
    Furthermore also the number length of the FIR filter is needed. The default value is 500. If the normalized doppler
    frequency is lower than 0.1 an equivalent process of higher normalized doppler is generated and linearly interpolated.

    References:

    [Stüber] Gordon L. Stüber, Principles of mobile communication, 2nd. ed., Kluwer, 2001.

    [Rappaport] Theodore S. Rappaport, Wireless communications: principles and practise, Prentice Hall, 1996.
  */
  class FIR_Fading_Generator : public Fading_Generator {
  public:
    // Set normalized dopper, doppler spectrum, length of FIR filter and calculation method
    FIR_Fading_Generator(const double norm_doppler = 0.0, const DOPPLER_SPECTRUM spectrum = Jakes, const int filter_length = 500);
    // Destructor
    virtual ~FIR_Fading_Generator(){ }

    //! Set FIR filter length
    void set_filter_length(const int filter_length);
    //! Get filter length
    int get_filter_length();

    //! Initialize the generator (is not needed)
    virtual void init();

    //! Generate \c no_samples values from the fading process
    virtual void generate(const int no_samples, cvec &output);
    //! Generate \c no_samples values from the fading process upsampled by \c upsampling_factor
    virtual void generate(const int no_samples, const int upsampling_factor, cvec &output);
    //! Generate \c no_samples values from the fading process
    //virtual cvec generate(const int no_samples);

  protected:
    int fir_length; // size of FIR filter
    int upsample_rate; // upsampling rate for linear interpolation
    MA_Filter<std::complex<double>, double, std::complex<double> > fir_filter;
    cvec left_overs; // left-overs from upsampling
  };

  /*!
    \ingroup channels
    \brief IFFT type Fading generator class
    \author Tony Ottosson

    A IFFT generator is a frequency domain filter implementation of filter method for generation of a fading process.
    Parameters that define the generator is the normalized doppler and the
    doppler spectrum. Possible values of doppler spectra are:
    <ul>
    <li> \a Jakes, J, Classic, or C: for the classical Jakes spectrum shape</li>
    <li> \a Rice or R: the classical Jakes spectrum and a direct tap. Observe that the LOS parameters also need to be set.</li>
    <li> \a GaussI, GI, Gauss1, or G1</li>
    <li> \a GaussII, GII, Gauss2, or G2</li>
    </ul>
    The method is block-based and consecutive blocks are independent fading. The method is very efficient for large blocks. 
    The size of the FFT, \f$N_\mathrm{fft}\f$, is given by the nearest higher integer power of two that is larger than \c no_samples.
    For small blocks the FFT size is increased to keep a good accuracy (at least 10 samples in the representation of the doppler-spectrum).
    However, to keep the program reasonably efficient the largest upsampling factor is 64. Higher factors will result in a run-time error.
    If so, please use another method.

    References:

    [Stüber] Gordon L. Stüber, Principles of mobile communication, 2nd. ed., Kluwer, 2001.

    [Rappaport] Theodore S. Rappaport, Wireless communications: principles and practise, Prentice Hall, 1996.
  */
  class IFFT_Fading_Generator : public Fading_Generator {
  public:
    // Set normalized dopper, doppler spectrum, length of FIR filter and calculation method
    IFFT_Fading_Generator(const double norm_doppler = 0.0, const DOPPLER_SPECTRUM spectrum = Jakes);
    // Destructor
    virtual ~IFFT_Fading_Generator(){ }

    //! Initialize the generator (is not needed)
    virtual void init();

    //! Generate \c no_samples values from the fading process
    virtual void generate(const int no_samples, cvec &output);
    //! Generate \c no_samples values from the fading process upsampled by \c upsampling_factor
    virtual void generate(const int no_samples, const int upsampling_factor, cvec &output);
    //! Generate \c no_samples values from the fading process
    //virtual cvec generate(const int no_samples);

  protected:
    //! Generator for Jakes spectrum
    void generate_Jakes(const int no_samples, cvec &output);

  };

  /*!
    \ingroup channels
    \brief General specification of a time-domain multipath channel
    \author Tony Ottosson

    This class does NOT generate any channel values. It is only used to specify a channel and to help resampling it
    to fit the sample time \c Ts of your need. To generate channel coefficients use the Tapped-Delay Line (TDL) class 
    TDL_Channel.

    A time invariant (or at least wide-sense stationary) channel have an impulse response that can be modeled as:
    \f[
    h(t) = \sum_{k=0}^{N_\mathrm{taps}-1} a_k \exp (-j \theta_k ) \delta(t-\tau_k) ,
    \f]
    where \f$N_taps\f$ is the number of channel taps, \f$a_k\f$ is the average amplitude at delay \f$\tau_k\f$, and \f$\theta_k\f$ is the channel
    phase of the \f$k\f$th channel tap. The average power profile,  and the delay profiles are defined as:
    \f[
    \mathbf{a} = [a_0, a_1, \ldots, a_{N_\mathrm{taps}-1}]
    \f]
    and
    \f[
    \mathbf{\tau} = [\tau_0, \tau_1, \ldots, \tau_{N_\mathrm{taps}-1}],
    \f]
    respectively. We assume without loss of generality that \f$\tau_0 = 0\f$ and
    \f$\tau_0 < \tau_1 < \ldots < \tau_{N_\mathrm{taps}-1}\f$.

    To initialize the class the following parameters should be defined:
    <ul>
    <li> \a avg_power_dB Average power profile, \f$\mathbf{a}\f$, given in dB</li>
    <li> \a delay_profile Delay profile, \f$\mathbf{\tau}\f$, given in seconds. Should be sorted (increasing delay) and first delay has to be 0</li>
    </ul>
    
    It is also possible to specify a predefined channel. The existing are:
    <ul>
    <li>ITU Channel models:</li>
    <ul>
    <li>\a ITU_Vehicular_A, 6-tap channel</li>
    <li>\a ITU_Vehicular_B, 6-tap channel</li>
    <li>\a ITU_Pedestrian_A, 4-tap channel</li>
    <li>\a ITU_Pedestrian_B, 6-tap channel</li>
    </ul>
    <li>COST 207 Channel models (see [Pätzold] pp. 259-266):
    <ul>
    <li>\a COST207_RA, Rural area, 4-tap channel</li>
    <li>\a COST207_RA6, Rural area, 6-tap channel</li>
    <li>\a COST207_TU, Typical urban, 6-tap channel</li>
    <li>\a COST207_TU6alt, Typical urban, alternative 6-tap channel</li>
    <li>\a COST207_TU12, Typical urban, 12-tap channel</li>
    <li>\a COST207_TU12alt, Typical urban, alternative 12-tap channel</li>
    <li>\a COST207_BU, Bad urban, 6-tap channel</li>
    <li>\a COST207_BU6alt, Bad urban, alternative 6-tap channel</li>
    <li>\a COST207_BU12, Bad urban, 12-tap channel</li>
    <li>\a COST207_BU12alt, Bad urban, alternative 12-tap channel</li>
    <li>\a COST207_HT, Hilly terrain, 6-tap channel</li>
    <li>\a COST207_HT6alt, Hilly terrain, alternative 6-tap channel</li>
    <li>\a COST207_HT12, Hilly terrain, 12-tap channel</li>
    <li>\a COST207_HT12alt, Hilly terrain, alternative 12-tap channel</li>
    </ul>
    <li>COST 259 Channel models (see [3GPP TR 25.943]):
    <ul>
    <li>\a COST259_TUx, Typical urban, 20-tap channel</li>
    <li>\a COST259_RAx, Rural ara, 10-tap channel</li>
    <li>\a COST259_HTx, Hilly terrain, 20-tap channel</li>
    </ul>
    </ul>

    Before assigning a channel specification to your TDL_Channel class the channel need to be
    discretized. This is done by calling discretize(Ts) where \c Ts is the sample time.

    References:

    [Pätzold] Matthias Pätzold, Mobile fading channels, Wiley, 2002.

    [3GPP TR 25.943] Technical Specification Group Radio Access Networs; Deployment aspects. Version 5.1.0 (2002-06).
  */
  class Channel_Specification {
  public:

    //! Initialize the average power profile in dB, and power delay profile in seconds
    Channel_Specification(const vec &avg_power_dB="0", const vec &delay_prof="0");
    //! Initialize with predetermined channel profile
    Channel_Specification(const CHANNEL_PROFILE profile);
    //! Destructor
    virtual ~Channel_Specification() { }

    //! Set both average power profile in dB and power delay profile in seconds
    void set_channel_profile(const vec &avg_power_dB="0", const vec &delay_prof="0");
    //! Set channel profile to a predetermined profile
    void set_channel_profile(const CHANNEL_PROFILE profile);

    //! Set doppler spectrum for each tap in the channel profile. If not set default is Jakes.
    void set_doppler_spectrum(DOPPLER_SPECTRUM *tap_spectrum);
    //! Set doppler spectrum for tap \c tap_number in the channel profile.
    void set_doppler_spectrum(const int tap_number, DOPPLER_SPECTRUM tap_spectrum);

    /*! Set LOS component for the first tap (zero delay).
      Only possible if Rice is chosen as doppler spectrum. Relative power (Rice factor) and doppler relative the maximum doppler frequency.
    */
    void set_LOS(const double relative_power, const double relative_doppler);

    //! Get both average power profile in dB and power delay profile in seconds
    void get_channel_profile(vec &avg_power_dB, vec &delay_prof);

    //! Return power profile in dB
    vec get_avg_power_dB();
    //! Return delay profile in seconds
    vec get_delay_prof();
    //! Get doppler spectrum for tap \c index
    DOPPLER_SPECTRUM get_doppler_spectrum(const int index);
    //! Get LOS relative power (Rice factor) on first tap (zero delay). Only if fist tap is of type Rice spectrum.
    double get_LOS_power();
    //! Get LOS doppler (relative to the maximum doppler) on first tap (zero delay). Only if fist tap is of type Rice spectrum.
    double get_LOS_doppler();
    
    //! Return the number of channel taps
    int taps() { return N_taps; }

    //! Calculate mean excess delay in samples
    double calc_mean_excess_delay();
    //! Calculate RMS delay spread in samples
    double calc_rms_delay_spread();

    //! Return true if channel profile is discretized. False otherwise
    bool is_discrete() { return discrete; }
    //! Get discrete Ts value
    double get_discrete_Ts() { return discrete_Ts; }

    //! Discretize the channel profile with resolution \c Ts. All taps within ((i-0.5)Ts,(i+0.5)Ts] will belong to the ith discrete tap.
    void discretize(const double Ts);

  protected:
    vec a_prof_dB, d_prof;
    Array<DOPPLER_SPECTRUM> tap_doppler_spectrum; // doppler spectrum for each tap
    int N_taps;
    double los_dopp; // doppler on LOS-component
    double los_power; // power on LOS-component
    bool discrete; // true if channel profile has been dicretized
    double discrete_Ts; // the sample time of discretization
  };



  /*!
    \ingroup channels
    \brief Tapped Delay Line (TDL) channel model
    \author Tony Ottosson

    A time invariant (or at least wide-sense stationary) channel have an impulse response that can be modeled as:
    \f[
    h(t) = \sum_{k=0}^{N_\mathrm{taps}-1} a_k \exp (-j \theta_k ) \delta(t-\tau_k) ,
    \f]
    where \f$N_taps\f$ is the number of channel taps, \f$a_k\f$ is the average amplitude at delay \f$\tau_k\f$, and \f$\theta_k\f$ is the channel
    phase of the \f$k\f$th channel tap. The average power profile,  and the delay profiles are defined as:
    \f[
    \mathbf{a} = [a_0, a_1, \ldots, a_{N_\mathrm{taps}-1}]
    \f]
    and
    \f[
    \mathbf{\tau} = [\tau_0, \tau_1, \ldots, \tau_{N_\mathrm{taps}-1}]
    \f],
    respectively. We assume without loss of generality that \f$\tau_0 = 0\f$ and
    \f$\tau_0 < \tau_1 < \ldots < \tau_{N_\mathrm{taps}-1}\f$.

    To initialize the channel profile the following parameters should be defined:
    <ul>
    <li>avg_power_dB Average power profile, \f$\mathbf{a}\f$, given in dB</li>
    <li>delay_profile Delay profile, \f$\mathbf{\tau}\f$, given in samples. Should be sorted (increasing delay) and first delay has to be 0</li>
    </ul>

    The correlation in time is decided by the normalized doppler frequency and the doppler spectrum. The normalized doppler frequency is calculated
    as \f$f_\mathrm{max} T_s\f$, where \f$f_\mathrm{max}\f$ is the maximum doppler frequency and \f$T_s\f$ the sample duration.
    Depending on the parameter \c norm_doppler the following models are used:
    <ul>
    <li>\a norm_doppler = 0 gives perfect interleaving. Channel coefficients are generated as independent complex Gaussian numbers</li>
    <li>\a norm_doppler > 0 gives correlated fading. The generation process depends on the selected fading generation model</li>
    </ul>


    Two main types of generation methods exist: the filter method and the Rice method. In the filter method the fading process
    is generated as a filtering of complex Gaussian process to achieve a given doppler spectrum. Currently there are two filter implementations
    <ul>
    <li>\a FIR. Finite Impulse Response (FIR) filter. Rather slow and inaccurate. Use with care. </li>
    <li>\a IFFT. Filtering in the frequency-domain. Block-based. Consecutive blocks are independent fading. Very efficient for large blocks.</li>
    </ul>
    The preferred method is the Rice method that approximate the fading process as a sum of sinusoids. Currently there is only one implementation,
    the Rice Method of Exact Doppler Spread (Rice_MEDS) which also is the default choice.

    To summarize, the currently supported fading generation models are:
    <ul>
    <li>\a FIR, filter FIR method.</li>
    <li>\a IFFT, filter IFFT method.</li>
    <li>\a Rice_MEDS, Rice MEDS method</li>
    </ul>


    Example: Simulation of WCDMA
    \code
    #include "itpp/itcomm.h"

    int main() {

    Channel_Specification channel_spec(COST259_RAx); // select the COST259 Rural Area model
    channel_spec.discretize(130.2e-9); // sample the channel at half a chip rate (=0.5/3.84e6)

    TDL_Channel channel(channel_spec); // initialize with defined channel profile
    channel.set_norm_doppler(0.01); // set the normalized doppler

    cvec transmitted_signal;
    //
    // Your code here generates the transmitted signal
    //

    Array<cvec> coeff;
    //The used channel values is returned in the 'coeff' which is an Array of cvec.
    cvec received_signal = channel(transmitted_signal, coeff);
    }
    \endcode

    References:

    [Pätzold] Matthias Pätzold, Mobile fading channels, Wiley, 2002.

    [Stüber] Gordon L. Stüber, Principles of mobile communication, 2nd. ed., Kluwer, 2001.

    [Rappaport] Theodore S. Rappaport, Wireless communications: principles and practise, Prentice Hall, 1996.
  */
  class TDL_Channel {
  public:
    //! Set  normalized doppler, average power profile in dB, and power delay profile in samples
    TDL_Channel(const double norm_doppler=0.0, const vec &avg_power_dB="0", const ivec &delay_prof="0");
    //! Set with defined Channel_Specification (need to be discretized)
    TDL_Channel(Channel_Specification &channel_spec);
    //! Destructor
    virtual ~TDL_Channel() { }

    //! Set both average power profile in dB and power delay profile in samples
    void set_channel_profile(const vec &avg_power_dB="0", const ivec &delay_prof="0");
    //! Set channel profile to uniform with \c no_taps taps
    void set_channel_profile_uniform(const int no_taps);
    //! Set channel profile to exponential with \c delay_spread
    void set_channel_profile_exponential(const double delay_spread);
    //! Set channel profile to a defined Channel_Specification (need to be discretized)
    void set_channel_profile(Channel_Specification &channel_spec);

    //! Set normalized doppler rate
    void set_norm_doppler(const double norm_doppler);
    //! Set doppler spectrum for each tap in the channel profile. If not set default is Jakes.
    void set_doppler_spectrum(DOPPLER_SPECTRUM *tap_spectrum);

    //! Set LOS component for the first tap (zero delay). Only possible if Rice is chosen as doppler spectrum. Relative power (Rice factor) and normalized doppler.
    void set_LOS(const double relative_power, const double norm_doppler);

    //! Set the fading generation method to \c method. Default is Rice_MEDS.
    void set_generation_method(const FADING_GENERATION_METHOD method=Rice_MEDS);

    //! Return the number of channel taps
    int taps() { return N_taps; }
    //! Get both average power profile in dB and power delay profile in samples
    void get_channel_profile(vec &avg_power_dB, ivec &delay_prof);
    //! Return power profile in dB
    vec get_avg_power_dB();
    //! Return delay profile in samples
    ivec get_delay_prof();
    //! Return normalized doppler rate
    double get_norm_doppler();
    //! Get LOS relative power (Rice factor) on first tap (zero delay). Only if fist tap is of type Rice spectrum.
    double get_LOS_power();
    //! Get LOS doppler (relative to the maximum doppler) on first tap (zero delay). Only if fist tap is of type Rice spectrum.
    double get_LOS_doppler();
    //! Return fading generation method
    FADING_GENERATION_METHOD get_generation_method();

    //! Calculate mean excess delay in samples
    double calc_mean_excess_delay();
    //! Calculate RMS delay spread in samples
    double calc_rms_delay_spread();

    //! Initialize all fading generators (needs to be called before generate() or filter())
    void init();

    //! Generate \c no_samples values of the channel
    virtual void generate(const int no_samples, Array<cvec> &channel_coeff);
    //! Generate \c no_samples values of the channel. Returns cmat with one tap per column
    virtual void generate(const int no_samples, cmat &channel_coeff);
    //! Generate \c no_samples values from the fading process upsampled by \c upsampling_factor
    virtual void generate(const int no_samples, const int upsampling_factor, Array<cvec> &channel_coeff);
    //! Generate \c no_samples values from the fading process upsampled by \c upsampling_factor
    virtual void generate(const int no_samples, const int upsampling_factor, cmat &channel_coeff);

    //! Shift generator time offset by a number of samples
    void shift_time_offset(const int no_samples);

    //! Filter the \c input with the known channel values \c channel_coeff (e.g. from the generate function)
    void filter_known_channel(const cvec &input, cvec &output, const Array<cvec> &channel_coeff);
    //! Filter the \c input with the known channel values \c channel_coeff (e.g. from the generate function)
    void filter_known_channel(const cvec &input, cvec &output, const cmat &channel_coeff);

    //! Generate channel coefficients and filter the \c input. Return output and channel coefficients
    void filter(const cvec &input, cvec &output, Array<cvec> &channel_coeff);
    //! Generate channel coefficients and filter the \c input. Return output and channel coefficients
    void filter(const cvec &input, cvec &output, cmat &channel_coeff);
    //! Generate channel coefficients and filter the \c input. Return output and channel coefficients
    cvec filter(const cvec &input, Array<cvec> &channel_coeff);
    //! Generate channel coefficients and filter the \c input. Return output and channel coefficients
    cvec filter(const cvec &input, cmat &channel_coeff);
    //! Generate channel coefficients and filter the \c input. Only return the output
    void filter(const cvec &input, cvec &output);
    //! Generate channel coefficients and filter the \c input. Only return output
    cvec filter(const cvec &input);

    //! Generate channel coefficients and filter the \c input. Return output and channel coefficients
    void operator()(const cvec &input, cvec &output, Array<cvec> &channel_coeff);
    //! Generate channel coefficients and filter the \c input. Return output and channel coefficients
    void operator()(const cvec &input, cvec &output, cmat &channel_coeff);
    //! Generate channel coefficients and filter the \c input. Return output and channel coefficients
    cvec operator()(const cvec &input, Array<cvec> &channel_coeff);
    //! Generate channel coefficients and filter the \c input. Return output and channel coefficients
    cvec operator()(const cvec &input, cmat &channel_coeff);
    //! Generate channel coefficients and filter the \c input. Only return output
    cvec operator()(const cvec &input);

    //! Calculate impulse-response on the supplied channel coefficients (produced by the generate() function)
    void calc_impulse_response(const Array<cvec> &channel_coeff, Array<cvec> &impulse_response);

    //! Calculate frequency-response on the supplied channel coefficients (produced by the generate() function)
    void calc_frequency_response(const Array<cvec> &channel_coeff, Array<cvec> &frequency_response, const int fft_size);
    //! Calculate frequency-response on the supplied channel coefficients (produced by the generate() function)
    void calc_frequency_response(const cmat &channel_coeff, cmat &frequency_response, const int fft_size);


  protected:
    vec a_prof; // average amplitude of each tap
    ivec d_prof; // delay in samples for each tap
    int N_taps; // number of taps
    double n_dopp; // normalized doppler
    double los_dopp; // doppler on LOS-component
    double los_power; // power on LOS-component
    Array<DOPPLER_SPECTRUM> tap_doppler_spectrum; // doppler spectrum for each tap
    Array<Fading_Generator *> fading_gen; // fading generators
    FADING_GENERATION_METHOD method; // generation method
    bool init_flag; // Channel ready to produce data
  };



  /*!
    \ingroup channels
    \brief A Binary Symetric Channel with crossover probability p.

    Input and output are of type \c bvec with 0 and 1.
    Example:
    \code
    #include "itpp/itcomm.h"

    int main() {
    //Initiate the BSC with cross-over probability 0.1
    BSC bsc(0.1);

    bvec transmitted_bits = randb(100);
    bvec received_bits = bsc(transmitted_bits);
    }
    \endcode
  */
  class BSC {
  public:
    //! Class constructor. Sets the error probability to p
    BSC(double in_p=0.0) : u(0.0, 1.0) { p = in_p; };
      //! Set crossover (bit error) probability
      void set_prob(double in_p) { p = in_p; };
      //! Get crossover (bit error) probability
      double get_prob() const { return p; };
      //! Feed \c input through the BSC channel
      bvec operator()(const bvec &input);
  private:
      Uniform_RNG u;
      double p;
  };



  /*!
    \ingroup channels
    \brief Ordinary AWGN Channel for cvec or vec inputs and outputs.

    For real signals, the input parameter (\c noisevar) should be set to
    \f$N_0/2\f$, where \f$N_0\f$ is the noise spectral density.
    However, in case of complex signals, the input parameter (\c noisevar) 
    represents the total noise variance of both real and imaginary parts,
    and thus is equal to \f$N_0\f$.

    Example:
    \code
    #include "itpp/itcomm.h"

    int main() {

    //Initiate the AWGN_Channel class
    double noisevar = 0.1;
    AWGN_Channel awgn_channel(noisevar);

    //Initiate a QPSK-modulator, and generate the transmitted signal
    QPSK qpsk;
    bvec transmitted_bits = randb(20);
    cvec transmitted_signal = qpsk.modulate_bits(transmitted_bits);

    //Usage of the member operator ()
    cvec received_signal = awgn_channel(transmitted_signal);

    //Demodulate the bits
    bvec received_bits = qpsk.demodulate_bits(received_signal);
    }
    \endcode
  */
  class AWGN_Channel {
  public:
    //! Class constructor. Sets the noise variance (for complex-valued channels the sum of real and imaginary parts)
    AWGN_Channel(double noisevar = 0.0) { sigma = std::sqrt(noisevar); }
    //! Set noise variance (for complex-valued channels the sum of real and imaginary parts)
    void set_noise(double noisevar) { sigma = std::sqrt(noisevar); }
    //! Get noise variance (for complex-valued channels the sum of real and imaginary parts)
    double get_noise() { return sqr(sigma); }
    //! Feed the complex input \c input through the complex-valued AWGN channel
    cvec operator()(const cvec &input);
    //! Feed the input \c through the real-valued AWGN channel
    vec operator()(const vec &input);
  protected:
    //! Standard deviation of the AWGN
    double sigma;
  };



} //namespace itpp

#endif // __channel_h
