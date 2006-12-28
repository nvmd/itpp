/*!
 * \file
 * \brief Implementation of classes for random number generators
 * \author Tony Ottosson
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

#include <ctime>
#ifdef _MSC_VER
#  include <process.h>
#  define getpid _getpid
#else
#  include <unistd.h>
#endif

#include <itpp/base/random.h>
#include <itpp/base/math/elem_math.h>


namespace itpp { 

  int Random_Generator::left = 0;
  const unsigned int Random_Generator::MAXINT = 0xffffffff;   // largest value from randInt()
  const unsigned int Random_Generator::MAGIC  = 0x9908b0dfU;  // magic constant

  unsigned int Random_Generator::state[624];
  unsigned int *Random_Generator::pNext;

  //! Variable used to ensure proper seed initialization 
  static bool __Random_Generator_seed_is_initialized = false; 

  ///////////////////////////////////////////////
  // Random_Generator
  ///////////////////////////////////////////////

  Random_Generator::Random_Generator() : lastSeed( 4357U )
  {
    if (!__Random_Generator_seed_is_initialized){
      reset();
      __Random_Generator_seed_is_initialized = true;
    }
  }

  void Random_Generator::randomize()
  {
    lastSeed = static_cast<unsigned long>(time(0)) * getpid(); // not good enough if randomize is used within a short time-interval
    reset();
  }

  void Random_Generator::get_state(ivec &out_state)
  {
    out_state.set_size(625, false);
    for (int i=0; i<624; i++)
      out_state(i) = state[i];

    out_state(624) = left; // the number of elements left in state before reload
    //saved_pNext = pNext;
  }

  void Random_Generator::set_state(ivec &new_state)
  {
    it_assert(new_state.size()==625, "Random_Generator::set_state(): Not a valid state vector");

    for (int i=0; i<624; i++)
      state[i] = new_state(i);

    left = new_state(624);
    pNext = &state[624-left];
  }

  // Set the seed of the Global Random Number Generator
  void RNG_reset(unsigned long seed)
  {
    Random_Generator RNG;
    RNG.reset(seed);
  }

  // Set the seed of the Global Random Number Generator to the same as last time
  void RNG_reset()
  {
    Random_Generator RNG;
    RNG.reset();
  }

  // Set a random seed for the Global Random Number Generator
  void RNG_randomize()
  {
    Random_Generator RNG;
    RNG.randomize();
  }

  // Save current full state of generator in memory
  void RNG_get_state(ivec &state)
  {
    Random_Generator RNG;
    RNG.get_state(state);
  }

  // Resume the state saved in memory
  void RNG_set_state(ivec &state)
  {
    Random_Generator RNG;
    RNG.set_state(state);
  }

  ///////////////////////////////////////////////
  // I_Uniform_RNG
  ///////////////////////////////////////////////

  I_Uniform_RNG::I_Uniform_RNG(int min, int max)
  {
    setup(min, max);
  }

  void I_Uniform_RNG::setup(int min, int max)
  {
    if (min <= max) {
      lo = min;
      hi = max;
    }
    else {
      lo = max;
      hi = min;
    }
  }

  void I_Uniform_RNG::get_setup(int &min, int &max) const
  {
    min = lo;
    max = hi;
  }

  ivec I_Uniform_RNG::operator()(int n)
  {
    ivec vv(n);

    for (int i=0; i<n; i++)
      vv(i) = sample();

    return vv;
  }

  imat I_Uniform_RNG::operator()(int h, int w)
  {
    imat mm(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	mm(i,j) = sample();

    return mm;
  }

  ///////////////////////////////////////////////
  // Uniform_RNG
  ///////////////////////////////////////////////

  Uniform_RNG::Uniform_RNG(double min, double max)
  {
    setup(min, max);
  }

  void Uniform_RNG::setup(double min, double max)
  {
    if (min <= max) {
      lo_bound = min;
      hi_bound = max;
    }
    else {
      lo_bound = max;
      hi_bound = min;
    }
  }

  void Uniform_RNG::get_setup(double &min, double &max) const
  {
    min = lo_bound;
    max = hi_bound;
  }

  ///////////////////////////////////////////////
  // Exp_RNG
  ///////////////////////////////////////////////

  Exponential_RNG::Exponential_RNG(double lambda)
  {
    setup(lambda);
  }

  vec Exponential_RNG::operator()(int n)
  {
    vec vv(n);

    for (int i=0; i<n; i++)
      vv(i) = sample();

    return vv;
  }

  mat Exponential_RNG::operator()(int h, int w)
  {
    mat mm(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	mm(i,j) = sample();

    return mm;
  }

  ///////////////////////////////////////////////
  // Normal_RNG
  ///////////////////////////////////////////////

  void Normal_RNG::get_setup(double &meanval, double &variance) const
  {
    meanval = mean;
    variance = sigma*sigma;
  }

  ///////////////////////////////////////////////
  // Laplace_RNG
  ///////////////////////////////////////////////

  Laplace_RNG::Laplace_RNG(double meanval, double variance)
  {
    setup(meanval, variance);
  }

  void Laplace_RNG::setup(double meanval, double variance)
  {
    mean = meanval;
    var = variance;
  }

  void Laplace_RNG::get_setup(double &meanval, double &variance) const
  {
    meanval = mean;
    variance = var;
  }



  vec Laplace_RNG::operator()(int n)
  {
    vec vv(n);

    for (int i=0; i<n; i++)
      vv(i) = sample();

    return vv;
  }

  mat Laplace_RNG::operator()(int h, int w)
  {
    mat mm(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	mm(i,j) = sample();

    return mm;
  }

  ///////////////////////////////////////////////
  // AR1_Normal_RNG
  ///////////////////////////////////////////////

  AR1_Normal_RNG::AR1_Normal_RNG(double meanval, double variance, double rho)
  {
    mean = meanval;
    var = variance;
    r = rho;
    mem = 0.0;
    factr = -2.0 * var * (1.0 - rho*rho);
    odd = true;
  }

  void AR1_Normal_RNG::setup(double meanval, double variance, double rho)
  {
    mean = meanval;
    var = variance;
    r = rho;
    factr = -2.0 * var * (1.0 - rho*rho);
    mem = 0.0;
    odd = true;
  }

  void AR1_Normal_RNG::get_setup(double &meanval, double &variance, double &rho) const
  {
    meanval = mean;
    variance = var;
    rho = r;
  }

  vec AR1_Normal_RNG::operator()(int n)
  {
    vec vv(n);

    for (int i=0; i<n; i++)
      vv(i) = sample();

    return vv;
  }

  mat AR1_Normal_RNG::operator()(int h, int w)
  {
    mat mm(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	mm(i,j) = sample();

    return mm;
  }

  void AR1_Normal_RNG::reset()
  {
    mem = 0.0;
  }

  ///////////////////////////////////////////////
  // Weibull_RNG
  ///////////////////////////////////////////////

  Weibull_RNG::Weibull_RNG(double lambda, double beta)
  {
    setup(lambda, beta);
  }

  void Weibull_RNG::setup(double lambda, double beta)
  {
    l=lambda;
    b=beta;
    mean = gamma(1.0 + 1.0 / b) / l;
    var = gamma(1.0+2.0/b)/(l*l) - mean;
  }


  vec Weibull_RNG::operator()(int n)
  {
    vec vv(n);

    for (int i=0; i<n; i++)
      vv(i) = sample();

    return vv;
  }

  mat Weibull_RNG::operator()(int h, int w)
  {
    mat mm(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	mm(i,j) = sample();

    return mm;
  }

  ///////////////////////////////////////////////
  // Rayleigh_RNG
  ///////////////////////////////////////////////

  Rayleigh_RNG::Rayleigh_RNG(double sigma)
  {
    setup(sigma);
  }

  vec Rayleigh_RNG::operator()(int n)
  {
    vec vv(n);

    for (int i=0; i<n; i++)
      vv(i) = sample();

    return vv;
  }

  mat Rayleigh_RNG::operator()(int h, int w)
  {
    mat mm(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	mm(i,j) = sample();

    return mm;
  }

  ///////////////////////////////////////////////
  // Rice_RNG
  ///////////////////////////////////////////////

  Rice_RNG::Rice_RNG(double lambda, double beta)
  {
    setup(lambda, beta);
  }

  vec Rice_RNG::operator()(int n)
  {
    vec vv(n);

    for (int i=0; i<n; i++)
      vv(i) = sample();

    return vv;
  }

  mat Rice_RNG::operator()(int h, int w)
  {
    mat mm(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	mm(i,j) = sample();

    return mm;
  }

} // namespace itpp
