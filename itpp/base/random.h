/*!
 * \file
 * \brief Definition of classes for random number generators
 * \author Tony Ottosson and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#ifndef RANDOM_H
#define RANDOM_H

#include <itpp/base/random_dsfmt.h>
#include <itpp/base/operators.h>


namespace itpp
{

//! \addtogroup randgen

/*!
 * \brief Base class for random (stochastic) sources.
 * \ingroup randgen
 *
 * Random_Generator is a typedef of DSFMT class specialization using 19937
 * generation period.
 *
 * \sa DSFMT
 */
typedef DSFMT_19937_RNG Random_Generator;

//! \addtogroup randgen
//!@{

//! Set the seed of the Global Random Number Generator
void RNG_reset(unsigned int seed);
//! Set the seed of the Global Random Number Generator to the same as last time
void RNG_reset();
//! Set a random seed for the Global Random Number Generator
void RNG_randomize();
//! Save the current state of the generator in a vector
void RNG_get_state(ivec &state);
//! Resume the state of the generator from previously saved vector
void RNG_set_state(const ivec &state);
//!@}

/*!
  \brief Bernoulli distribution
  \ingroup randgen
*/
class Bernoulli_RNG
{
public:
  //! Binary source with probability prob for a 1
  Bernoulli_RNG(double prob) { setup(prob); }
  //! Binary source with probability prob for a 1
  Bernoulli_RNG() { p = 0.5; }
  //! set the probability
  void setup(double prob) {
    it_assert(prob >= 0.0 && prob <= 1.0, "The Bernoulli source probability "
              "must be between 0 and 1");
    p = prob;
  }
  //! return the probability
  double get_setup() const { return p; }
  //! Get one sample.
  bin operator()() { return sample(); }
  //! Get a sample vector.
  bvec operator()(int n) { bvec temp(n); sample_vector(n, temp); return temp; }
  //! Get a sample matrix.
  bmat operator()(int h, int w) { bmat temp(h, w); sample_matrix(h, w, temp); return temp; }
  //! Get a sample
  bin sample() { return RNG.genrand_close_open() < p ? bin(1) : bin(0); }
  //! Get a sample vector.
  void sample_vector(int size, bvec &out) {
    out.set_size(size, false);
    for (int i = 0; i < size; i++) out(i) = sample();
  }
  //! Get a sample matrix.
  void sample_matrix(int rows, int cols, bmat &out) {
    out.set_size(rows, cols, false);
    for (int i = 0; i < rows*cols; i++) out(i) = sample();
  }
protected:
private:
  //!
  double p;
  //!
  Random_Generator RNG;
};

/*!
  \brief Integer uniform distribution
  \ingroup randgen

  Example: Generation of random uniformly distributed integers in the interval [0,10].
  \code
  #include "itpp/sigproc.h"

  int main() {

  I_Uniform_RNG gen(0, 10);

  cout << gen() << endl; // prints a random integer
  cout << gen(10) << endl; // prints 10 random integers
  }
  \endcode
*/
class I_Uniform_RNG
{
public:
  //! constructor. Sets min and max values.
  I_Uniform_RNG(int min = 0, int max = 1);
  //! set min and max values
  void setup(int min, int max);
  //! get the parameters
  void get_setup(int &min, int &max) const;
  //! Get one sample.
  int operator()() { return sample(); }
  //! Get a sample vector.
  ivec operator()(int n);
  //! Get a sample matrix.
  imat operator()(int h, int w);
  //! Return a single value from this random generator
  int sample() {
    return floor_i(RNG.genrand_close_open() * (hi - lo + 1)) + lo;
  }
private:
  //!
  int lo;
  //!
  int hi;
  //!
  Random_Generator RNG;
};

/*!
  \brief Uniform distribution
  \ingroup randgen
*/
class Uniform_RNG
{
public:
  //! Constructor. Set min, max and seed.
  Uniform_RNG(double min = 0, double max = 1.0);
  //! set min and max
  void setup(double min, double max);
  //! get parameters
  void get_setup(double &min, double &max) const;
  //! Get one sample.
  double operator()() { return (sample() * (hi_bound - lo_bound) + lo_bound); }
  //! Get a sample vector.
  vec operator()(int n) {
    vec temp(n);
    sample_vector(n, temp);
    temp *= hi_bound - lo_bound;
    temp += lo_bound;
    return temp;
  }
  //! Get a sample matrix.
  mat operator()(int h, int w) {
    mat temp(h, w);
    sample_matrix(h, w, temp);
    temp *= hi_bound - lo_bound;
    temp += lo_bound;
    return temp;
  }
  //! Get a Uniformly distributed [0,1) sample
  double sample() {  return RNG.genrand_close_open(); }
  //! Get a Uniformly distributed [0,1) vector
  void sample_vector(int size, vec &out) {
    out.set_size(size, false);
    for (int i = 0; i < size; i++) out(i) = sample();
  }
  //! Get a Uniformly distributed [0,1) matrix
  void sample_matrix(int rows, int cols, mat &out) {
    out.set_size(rows, cols, false);
    for (int i = 0; i < rows*cols; i++) out(i) = sample();
  }
protected:
private:
  //!
  double lo_bound, hi_bound;
  //!
  Random_Generator RNG;
};

/*!
  \brief Exponential distribution
  \ingroup randgen
*/
class Exponential_RNG
{
public:
  //! constructor. Set lambda.
  Exponential_RNG(double lambda = 1.0);
  //! Set lambda
  void setup(double lambda) { l = lambda; }
  //! get lambda
  double get_setup() const;
  //! Get one sample.
  double operator()() { return sample(); }
  //! Get a sample vector.
  vec operator()(int n);
  //! Get a sample matrix.
  mat operator()(int h, int w);
private:
  double sample() { return (-std::log(RNG.genrand_open_close()) / l); }
  double l;
  Random_Generator RNG;
};

/*!
 * \brief Normal distribution
 * \ingroup randgen
 *
 * Normal (Gaussian) random variables, using a simplified Ziggurat method.
 *
 * For details see the following arcticle: George Marsaglia, Wai Wan
 * Tsang, "The Ziggurat Method for Generating Random Variables", Journal
 * of Statistical Software, vol. 5 (2000), no. 8
 *
 * This implementation is based on the generator written by Jochen Voss
 * found at http://seehuhn.de/comp/ziggurat/, which is also included in
 * the GSL library (randlist/gauss.c).
 */
class Normal_RNG
{
public:
  //! Constructor. Set mean and variance.
  Normal_RNG(double meanval, double variance):
      mean(meanval), sigma(std::sqrt(variance)) {}
  //! Constructor. Set mean and variance.
  Normal_RNG(): mean(0.0), sigma(1.0) {}
  //! Set mean, and variance
  void setup(double meanval, double variance)
  { mean = meanval; sigma = std::sqrt(variance); }
  //! Get mean and variance
  void get_setup(double &meanval, double &variance) const;
  //! Get one sample.
  double operator()() { return (sigma*sample() + mean); }
  //! Get a sample vector.
  vec operator()(int n) {
    vec temp(n);
    sample_vector(n, temp);
    temp *= sigma;
    temp += mean;
    return temp;
  }
  //! Get a sample matrix.
  mat operator()(int h, int w) {
    mat temp(h, w);
    sample_matrix(h, w, temp);
    temp *= sigma;
    temp += mean;
    return temp;
  }
  //! Get a Normal distributed (0,1) sample
  double sample();

  //! Get a Normal distributed (0,1) vector
  void sample_vector(int size, vec &out) {
    out.set_size(size, false);
    for (int i = 0; i < size; i++) out(i) = sample();
  }

  //! Get a Normal distributed (0,1) matrix
  void sample_matrix(int rows, int cols, mat &out) {
    out.set_size(rows, cols, false);
    for (int i = 0; i < rows*cols; i++) out(i) = sample();
  }
private:
  double mean, sigma;
  static const double ytab[128];
  static const unsigned int ktab[128];
  static const double wtab[128];
  static const double PARAM_R;
  Random_Generator RNG;
};

/*!
 * \brief Gamma distribution
 * \ingroup randgen
 *
 * Generate samples from Gamma(alpha,beta) density, according to the
 * following equation:
 * \f[ x \sim \Gamma(\alpha,\beta) =
 * \frac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha-1} \exp(-\beta x) \f]
 *
 * For \f$\alpha=1\f$ the Gamma distribution is equivalent to the
 * Exponential distribution.
 *
 * \note The implementation of the sample() function was adapted from the
 * R statistical language.
 * \author Vasek Smidl
 */
class Gamma_RNG
{
public:
  //! Constructor, which sets alpha (a) and beta (b)
  Gamma_RNG(double a = 1.0, double b = 1.0): alpha(a), beta(b) {}
  //! Set alpha and beta
  void setup(double a, double b) { alpha = a; beta = b; }
  //! Get one sample
  double operator()() { return sample(); }
  //! Get a sample vector
  vec operator()(int n);
  //! Get a sample matrix
  mat operator()(int r, int c);
  //! Get a sample
  double sample();
private:
  double alpha;
  double beta;
  Random_Generator RNG;
  Normal_RNG NRNG;
};

/*!
  \brief Laplacian distribution
  \ingroup randgen
*/
class Laplace_RNG
{
public:
  //! Constructor. Set mean and variance.
  Laplace_RNG(double meanval = 0.0, double variance = 1.0);
  //! Set mean and variance
  void setup(double meanval, double variance);
  //! Get mean and variance
  void get_setup(double &meanval, double &variance) const;
  //! Get one sample.
  double operator()() { return sample(); }
  //! Get a sample vector.
  vec operator()(int n);
  //! Get a sample matrix.
  mat operator()(int h, int w);
  //! Returns a single sample
  double sample() {
    double u = RNG.genrand_open_open();
    double l = sqrt_12var;
    if (u < 0.5)
      l *= std::log(2.0 * u);
    else
      l *= -std::log(2.0 * (1 - u));
    return (l + mean);
  }
private:
  double mean, var, sqrt_12var;
  Random_Generator RNG;
};

/*!
  \brief A Complex Normal Source
  \ingroup randgen
*/
class Complex_Normal_RNG
{
public:
  //! Constructor. Set mean and variance.
  Complex_Normal_RNG(std::complex<double> mean, double variance):
      norm_factor(1.0 / std::sqrt(2.0)) {
    setup(mean, variance);
  }
  //! Default constructor
  Complex_Normal_RNG(): m(0.0), sigma(1.0), norm_factor(1.0 / std::sqrt(2.0)) {}
  //! Set mean and variance
  void setup(std::complex<double> mean, double variance) {
    m = mean;
    sigma = std::sqrt(variance);
  }
  //! Get mean and variance
  void get_setup(std::complex<double> &mean, double &variance) {
    mean = m;
    variance = sigma * sigma;
  }
  //! Get one sample.
  std::complex<double> operator()() { return sigma*sample() + m; }
  //! Get a sample vector.
  cvec operator()(int n) {
    cvec temp(n);
    sample_vector(n, temp);
    temp *= sigma;
    temp += m;
    return temp;
  }
  //! Get a sample matrix.
  cmat operator()(int h, int w) {
    cmat temp(h, w);
    sample_matrix(h, w, temp);
    temp *= sigma;
    temp += m;
    return temp;
  }
  //! Get a Complex Normal (0,1) distributed sample
  std::complex<double> sample() {
    double a = nRNG.sample() * norm_factor;
    double b = nRNG.sample() * norm_factor;
    return std::complex<double>(a, b);
  }

  //! Get a Complex Normal (0,1) distributed vector
  void sample_vector(int size, cvec &out) {
    out.set_size(size, false);
    for (int i = 0; i < size; i++) out(i) = sample();
  }

  //! Get a Complex Normal (0,1) distributed matrix
  void sample_matrix(int rows, int cols, cmat &out) {
    out.set_size(rows, cols, false);
    for (int i = 0; i < rows*cols; i++) out(i) = sample();
  }

  //! Dummy assignment operator - MSVC++ warning C4512
  Complex_Normal_RNG & operator=(const Complex_Normal_RNG&) { return *this; }

private:
  std::complex<double> m;
  double sigma;
  const double norm_factor;
  Normal_RNG nRNG;
};

/*!
  \brief Filtered normal distribution
  \ingroup randgen
*/
class AR1_Normal_RNG
{
public:
  //! Constructor. Set mean, variance, and correlation.
  AR1_Normal_RNG(double meanval = 0.0, double variance = 1.0,
                 double rho = 0.0);
  //! Set mean, variance, and correlation
  void setup(double meanval, double variance, double rho);
  //! Get mean, variance and correlation
  void get_setup(double &meanval, double &variance, double &rho) const;
  //! Set memory contents to zero
  void reset();
  //! Get a single random sample
  double operator()() { return sample(); }
  //! Get a sample vector.
  vec operator()(int n);
  //! Get a sample matrix.
  mat operator()(int h, int w);
private:
  double sample() {
    mem *= r;
    if (odd) {
      r1 = m_2pi * RNG.genrand_open_close();
      r2 = std::sqrt(factr * std::log(RNG.genrand_open_close()));
      mem += r2 * std::cos(r1);
    }
    else {
      mem += r2 * std::sin(r1);
    }
    odd = !odd;
    return (mem + mean);
  }
  double mem, r, factr, mean, var, r1, r2;
  bool odd;
  Random_Generator RNG;
};

/*!
  \brief Gauss_RNG is the same as Normal Source
  \ingroup randgen
*/
typedef Normal_RNG Gauss_RNG;

/*!
  \brief AR1_Gauss_RNG is the same as AR1_Normal_RNG
  \ingroup randgen
*/
typedef AR1_Normal_RNG AR1_Gauss_RNG;

/*!
  \brief Weibull distribution
  \ingroup randgen
*/
class Weibull_RNG
{
public:
  //! Constructor. Set lambda and beta.
  Weibull_RNG(double lambda = 1.0, double beta = 1.0);
  //! Set lambda, and beta
  void setup(double lambda, double beta);
  //! Get lambda and beta
  void get_setup(double &lambda, double &beta) { lambda = l; beta = b; }
  //! Get one sample.
  double operator()() { return sample(); }
  //! Get a sample vector.
  vec operator()(int n);
  //! Get a sample matrix.
  mat operator()(int h, int w);
private:
  double sample() {
    return (std::pow(-std::log(RNG.genrand_open_close()), 1.0 / b) / l);
  }
  double l, b;
  double mean, var;
  Random_Generator RNG;
};

/*!
  \brief Rayleigh distribution
  \ingroup randgen
*/
class Rayleigh_RNG
{
public:
  //! Constructor. Set sigma.
  Rayleigh_RNG(double sigma = 1.0);
  //! Set sigma
  void setup(double sigma) { sig = sigma; }
  //! Get sigma
  double get_setup() { return sig; }
  //! Get one sample.
  double operator()() { return sample(); }
  //! Get a sample vector.
  vec operator()(int n);
  //! Get a sample matrix.
  mat operator()(int h, int w);
private:
  double sample() {
    double s1 = nRNG.sample();
    double s2 = nRNG.sample();
    // s1 and s2 are N(0,1) and independent
    return (sig * std::sqrt(s1*s1 + s2*s2));
  }
  double sig;
  Normal_RNG nRNG;
};

/*!
  \brief Rice distribution
  \ingroup randgen
*/
class Rice_RNG
{
public:
  //! Constructor. Set sigma, and v (if v = 0, Rice -> Rayleigh).
  Rice_RNG(double sigma = 1.0, double v = 1.0);
  //! Set sigma, and v (if v = 0, Rice -> Rayleigh).
  void setup(double sigma, double v) { sig = sigma; s = v; }
  //! Get parameters
  void get_setup(double &sigma, double &v) { sigma = sig; v = s; }
  //! Get one sample
  double operator()() { return sample(); }
  //! Get a sample vector
  vec operator()(int n);
  //! Get a sample matrix
  mat operator()(int h, int w);
private:
  double sample() {
    double s1 = nRNG.sample() + s;
    double s2 = nRNG.sample();
    // s1 and s2 are N(0,1) and independent
    return (sig * std::sqrt(s1*s1 + s2*s2));
  }
  double sig, s;
  Normal_RNG nRNG;
};

//! \addtogroup randgen
//!@{

//! Generates a random bit (equally likely 0s and 1s)
inline bin randb(void) { Bernoulli_RNG src; return src.sample(); }
//! Generates a random bit vector (equally likely 0s and 1s)
inline void randb(int size, bvec &out) { Bernoulli_RNG src; src.sample_vector(size, out); }
//! Generates a random bit vector (equally likely 0s and 1s)
inline bvec randb(int size) { bvec temp; randb(size, temp); return temp; }
//! Generates a random bit matrix (equally likely 0s and 1s)
inline void randb(int rows, int cols, bmat &out) { Bernoulli_RNG src; src.sample_matrix(rows, cols, out); }
//! Generates a random bit matrix (equally likely 0s and 1s)
inline bmat randb(int rows, int cols) { bmat temp; randb(rows, cols, temp); return temp; }

//! Generates a random uniform (0,1) number
inline double randu(void) { Uniform_RNG src; return src.sample(); }
//! Generates a random uniform (0,1) vector
inline void randu(int size, vec &out) { Uniform_RNG src; src.sample_vector(size, out); }
//! Generates a random uniform (0,1) vector
inline vec randu(int size) { vec temp; randu(size, temp); return temp; }
//! Generates a random uniform (0,1) matrix
inline void randu(int rows, int cols, mat &out) { Uniform_RNG src; src.sample_matrix(rows, cols, out); }
//! Generates a random uniform (0,1) matrix
inline mat randu(int rows, int cols) { mat temp; randu(rows, cols, temp); return temp; }

//! Generates a random integer in the interval [low,high]
inline int randi(int low, int high) { I_Uniform_RNG src; src.setup(low, high); return src(); }
//! Generates a random ivec with elements in the interval [low,high]
inline ivec randi(int size, int low, int high) { I_Uniform_RNG src; src.setup(low, high); return src(size); }
//! Generates a random imat with elements in the interval [low,high]
inline imat randi(int rows, int cols, int low, int high) { I_Uniform_RNG src; src.setup(low, high); return src(rows, cols); }

//! Generates a random Rayleigh vector
inline vec randray(int size, double sigma = 1.0) { Rayleigh_RNG src; src.setup(sigma); return src(size); }

//! Generates a random Rice vector (See J.G. Poakis, "Digital Communications, 3rd ed." p.47)
inline vec randrice(int size, double sigma = 1.0, double s = 1.0) { Rice_RNG src; src.setup(sigma, s); return src(size); }

//! Generates a random complex Gaussian vector
inline vec randexp(int size, double lambda = 1.0) { Exponential_RNG src; src.setup(lambda); return src(size); }

//! Generates a random Gaussian (0,1) variable
inline double randn(void) { Normal_RNG src; return src.sample(); }
//! Generates a random Gaussian (0,1) vector
inline void randn(int size, vec &out) { Normal_RNG src; src.sample_vector(size, out); }
//! Generates a random Gaussian (0,1) vector
inline vec randn(int size) { vec temp; randn(size, temp); return temp; }
//! Generates a random Gaussian (0,1) matrix
inline void randn(int rows, int cols, mat &out)  { Normal_RNG src; src.sample_matrix(rows, cols, out); }
//! Generates a random Gaussian (0,1) matrix
inline mat randn(int rows, int cols) { mat temp; randn(rows, cols, temp); return temp; }

/*! \brief Generates a random complex Gaussian (0,1) variable

The real and imaginary parts are independent and have variances equal to 0.5
*/
inline std::complex<double> randn_c(void) { Complex_Normal_RNG src; return src.sample(); }
//! Generates a random complex Gaussian (0,1) vector
inline void randn_c(int size, cvec &out)  { Complex_Normal_RNG src; src.sample_vector(size, out); }
//! Generates a random complex Gaussian (0,1) vector
inline cvec randn_c(int size) { cvec temp; randn_c(size, temp); return temp; }
//! Generates a random complex Gaussian (0,1) matrix
inline void randn_c(int rows, int cols, cmat &out) { Complex_Normal_RNG src; src.sample_matrix(rows, cols, out); }
//! Generates a random complex Gaussian (0,1) matrix
inline cmat randn_c(int rows, int cols) { cmat temp; randn_c(rows, cols, temp); return temp; }

//!@}

} // namespace itpp

#endif // #ifndef RANDOM_H
