/*!
 * \file
 * \brief Definitions of statistics functions and classes
 * \author Tony Ottosson, Johan Bergman, Adam Piatyszeka and Andy Panov
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

#ifndef STAT_H
#define STAT_H

#include <cmath>
#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/elmatfunc.h>
#include <itpp/base/matfunc.h>


namespace itpp {

  //! \addtogroup statistics
  //!@{

  /*! 
    \brief A class for sampling a signal and calculating statistics
  */
  class Stat {
  public:
    //! Default constructor
    Stat() {clear();}
    //! Destructor
    virtual ~Stat() {}

    //! Clear statistics
    virtual void clear()
    {
      _n_overflows = 0;
      _n_samples = 0;
      _n_zeros = 0;
      _max = 0.0;
      _min = 0.0;
      _sqr_sum = 0.0;
      _sum = 0.0;
    }

    //! Register a sample and flag for overflow
    virtual void sample(const double s, const bool overflow=false)
    {
      _n_samples++;
      _sum += s;
      _sqr_sum += s*s;
      if (s < _min) _min = s;
      if (s > _max) _max = s;
      if (overflow) _n_overflows++;
      if (s == 0.0) _n_zeros++;
    }

    //! Number of reported overflows
    int n_overflows() const {return _n_overflows;}
    //! Number of samples
    int n_samples() const {return _n_samples;}
    //! Number of zero samples
    int n_zeros() const {return _n_zeros;}
    //! Average over all samples
    double avg() const {return _sum/_n_samples;}
    //! Maximum sample
    double max() const {return _max;}
    //! Minimum sample
    double min() const {return _min;}
    //! Standard deviation of all samples
    double sigma() const
    {
      double sigma2 = _sqr_sum/_n_samples - avg()*avg();
      return std::sqrt(sigma2 < 0 ? 0 : sigma2);
    }
    //! Squared sum of all samples
    double sqr_sum() const {return _sqr_sum;}
    //! Sum of all samples
    double sum() const {return _sum;}
    //! Histogram over all samples (not implemented yet)
    vec histogram() const {return vec(0);}

  protected:
    //! Number of reported overflows
    int _n_overflows;
    //! Number of samples
    int _n_samples;
    //! Number of zero samples
    int _n_zeros;
    //! Maximum sample
    double _max;
    //! Minimum sample
    double _min;
    //! Squared sum of all samples
    double _sqr_sum;
    //! Sum of all samples
    double _sum;
  };


  // ------------------------------------------------------------------------

  /*! 
    \brief Histogram computation class
    \author Andy Panov
   
    The Histogram class counts the number of observations of arbitrary
    numerical types that fall into specified bins. Centers of the leftmost
    and rightmost bin along with a total number of bins are passed to a
    histogram object as constructor parameters. Histogram counters are
    updated when calling update() method. It is possible to access bin
    counters and bin interval parameters for all bins at once or separately
    for each bin.
      
    Example:
    \code
    // Create histogram with 100 bins spanning from 0 to 99 (leftmost bin is 
    // centered at 0, rightmost bin is centerd at 99). 
    Histogram<double> hist(0, 99, 100);
      
    // Compute histogram of 100 random variables taken from normal distribution
    hist.update(randn(100));
   
    // Get position of bin number 5
    double bin5_center = hist.get_bin_center(5);
    // Get corresponding bin counter
    int bin5_counter = hist.get_bin(5);
    // Get bin 5 left boundary:
    double bin5_left = hist.get_bin_left(5);
   
    // compute PDF & CDF of experimental data
    vec my_data_pdf = hist.get_pdf();
    vec my_data_cdf = hist.get_cdf();
    \endcode
   */
  template<typename Num_T>
  class Histogram {
  public:
    //! Default constructor. Constructs histogram with 100 bins spanning
    //! values from 0 to 99 by default.
    Histogram(Num_T from = Num_T(0), Num_T to = Num_T(99), int n_bins = 100);
    //! Default destructor
    ~Histogram() {};

    //! Histogram setup
    void setup(Num_T from, Num_T to, int n_bins);
    
    //! Histogram update
    void update(Num_T value);
    void update(Vec<Num_T> values);
    void update(Mat<Num_T> values);
    
    //! Bins reset, so accumulation can be restarted
    void reset() { trials_cnt = 0; bins.zeros(); };
    //! Access to single bin counter
    int get_bin(int ix) const { return bins(ix); };
    //! Access to histogram as a vector
    ivec get_bins() const { return bins; };
    //! Access to bin center values (all bins)
    Vec<Num_T> get_bin_centers() const { return center_vals; };
    //! Access to bin center (single bin)
    Num_T get_bin_center(int ix) const { return center_vals(ix); };
    //! Access to left boundary of bin intervals (all bins)
    Vec<Num_T> get_bin_lefts() const { return lo_vals; };
    //! Access to left boundary of single bin 
    Num_T get_bin_left(int ix) const { return lo_vals(ix); };
    //! Access to right boundary of bin intervals (all bins)
    Vec<Num_T> get_bin_rights() const { return hi_vals; };
    //! Access to right boundary of single bin 
    Num_T get_bin_right(int ix) const { return hi_vals(ix); };
    
    //! Experimental Probability Density Function (PDF) computation
    vec get_pdf() const;
    //! Experimental Cumulative Density Function (CDF) computation	
    vec get_cdf() const;
    
    //! Current number of bins
    int bins_num() const { return num_bins; };
    //! Current trials counter
    int trials_num() const {return trials_cnt;};
  
  private:
    //! Number of bins
    int num_bins;
    //! Step between bins
    Num_T step;
    //! Low boundaries of histogram bins
    Vec<Num_T> lo_vals;
    //! Upper boundaries of histogram bins
    Vec<Num_T> hi_vals;
    //! Bin centers
    Vec<Num_T> center_vals; 
    //! Bins storage
    ivec bins;
    //! Number of processed samples
    int trials_cnt;
  };

  template<class Num_T>
  inline Histogram<Num_T>::Histogram(Num_T from, Num_T to, int n_bins)
    
  {
    setup(from, to, n_bins);
  }

  template<class Num_T>
  inline void Histogram<Num_T>::setup(Num_T from, Num_T to, int n_bins)
  {
    num_bins = n_bins;
    lo_vals.set_size(n_bins);
    hi_vals.set_size(n_bins);
    center_vals.set_size(n_bins);
    bins.set_size(n_bins);
    trials_cnt = 0;
    step = (to - from) / (num_bins - 1);
    center_vals = linspace(from, to, num_bins);
    lo_vals = center_vals - step/2; 
    hi_vals = center_vals + step/2; 
    reset();
  }
	
  template<class Num_T>
  inline void Histogram<Num_T>::update(Num_T value) 
  {		
    // search for the corresponding bin using dichotomy approach
    int start = 0;
    int end = num_bins - 1;
    int test = (start + end) / 2;

    while (start < end) {
      if (value < lo_vals(test))
	end = test - 1;
      else if (value >= hi_vals(test))
	start = test + 1;
      else
	break;
      test = (start + end) / 2;
    };

    bins(test)++;
    trials_cnt++;
  }

 template<class Num_T>
 inline void Histogram<Num_T>::update(Vec<Num_T> values) 
 {		
   for (int i = 0; i < values.length(); i++)
     update(values(i));
 }

 template<class Num_T>
 inline void Histogram<Num_T>::update(Mat<Num_T> values) 
 {		
   for (int i = 0; i < values.rows(); i++)
     for (int j = 0; j < values.cols(); j++)
       update(values(i,j));            
 }

  template<class Num_T>
  inline vec Histogram<Num_T>::get_pdf() const 
  {
    vec pdf(num_bins);
    for (int j = 0; j < num_bins; j++)
      pdf(j) = static_cast<double>(bins(j)) / trials_cnt;
    return pdf;
  }

  template<class Num_T>
  inline vec Histogram<Num_T>::get_cdf() const 
  {
    ivec tmp = cumsum(bins);			
    vec cdf(num_bins);
    for (int j = 0; j < num_bins; j++)
      cdf(j) = static_cast<double>(tmp(j)) / trials_cnt;		
    return cdf;
  }

  // ------------------------------------------------------------------------

  
  //! Maximum value of vector
  template<class T>
    T max(const Vec<T> &v)
  {
    T	maxdata = v(0);
    for (int i=1; i<v.length(); i++)
      if (v(i)>maxdata)
	maxdata=v(i);
    return maxdata;
  }

  //! Maximum value of vector, also returns the index position of max value
  template<class T>
    T max(const Vec<T> &v, int& index)
  {
    T maxdata = v(0);
    index = 0;
    for (int i=1; i<v.length(); i++)
      if (v(i)>maxdata) {
	maxdata=v(i);
	index = i;
      }
    return maxdata;
  }

  /*! 
   * Maximum values over each row/column in the matrix \c m
   *
   * <tt>max(m) = max(m, 1)</tt> returns a vector where the elements are
   * maximum over each column, whereas <tt>max(m, 2)</tt> returns a vector
   * where the elements are maximum over each row.
   */
  template<class T>
    Vec<T> max(const Mat<T> &m, int dim=1)
  {
    it_assert(dim==1 || dim==2, "max: dimension need to be 1 or 2");
    Vec<T> out;

    if (dim == 1) {
      out.set_size(m.cols(), false);

      for (int i=0; i<m.cols(); i++)
	out(i) = max(m.get_col(i));
    } else {
      out.set_size(m.rows(), false);

      for (int i=0; i<m.rows(); i++)
	out(i) = max(m.get_row(i));
    }
      
    return out;
  }

  /*! 
   * Maximum values over each row/column in the matrix \c m
   *
   * <tt>max(m) = max(m, 1)</tt> returns a vector where the elements are
   * maximum over each column, whereas <tt>max(m, 2)</tt> returns a vector
   * where the elements are maximum over each row.
   *
   * Also returns a vector of indices with positions of maximum value within
   * a column/row.
   */
  template<class T>
    Vec<T> max(const Mat<T> &m, ivec &index, int dim=1)
  {
    it_assert(dim==1 || dim==2, "max: dimension need to be 1 or 2");
    Vec<T> out;

    if (dim == 1) {
      out.set_size(m.cols(), false);
      index.set_size(m.cols(), false);

      for (int i=0; i<m.cols(); i++)
	out(i) = max(m.get_col(i), index(i));
    } else {
      out.set_size(m.rows(), false);
      index.set_size(m.rows(), false);

      for (int i=0; i<m.rows(); i++)
	out(i) = max(m.get_row(i), index(i));
    }
      
    return out;
  }

  //! Minimum value of vector
  template<class T>
    T min(const Vec<T> &in)
  {
    T mindata=in[0];
    for (int i=1;i<in.length();i++)
      if (in[i]<mindata)
	mindata=in[i];
    return mindata;
  }

  //! Minimum value of vector, also returns the index position of min value
  template<class T>
    T min(const Vec<T> &in, int& index)
  {
    T mindata=in[0];
    index = 0;
    for (int i=1;i<in.length();i++)
      if (in[i]<mindata) {
	mindata=in[i];
	index = i;
      }
    return mindata;
  }


  /*! 
   * Minimum values over each row/column in the matrix \c m
   *
   * <tt>min(m) = min(m, 1)</tt> returns a vector where the elements are
   * minimum over each column, whereas <tt>min(m, 2)</tt> returns a vector
   * where the elements are minimum over each row.
   */
  template<class T>
    Vec<T> min(const Mat<T> &m, int dim=1)
  {
    it_assert(dim==1 || dim==2, "min: dimension need to be 1 or 2");
    Vec<T> out;

    if (dim == 1) {
      out.set_size(m.cols(), false);

      for (int i=0; i<m.cols(); i++)
	out(i) = min(m.get_col(i));
    } else {
      out.set_size(m.rows(), false);

      for (int i=0; i<m.rows(); i++)
	out(i) = min(m.get_row(i));
    }
    return out;
  }


  /*! 
   * Minimum values over each row/column in the matrix \c m
   *
   * <tt>min(m) = min(m, 1)</tt> returns a vector where the elements are
   * minimum over each column, whereas <tt>min(m, 2)</tt> returns a vector
   * where the elements are minimum over each row.
   *
   * Also returns a vector of indices with positions of minimum value within
   * a column/row.
   */
  template<class T>
    Vec<T> min(const Mat<T> &m,  ivec &index, int dim=1)
  {
    it_assert(dim==1 || dim==2, "min: dimension need to be 1 or 2");
    Vec<T> out;

    if (dim == 1) {
      out.set_size(m.cols(), false);
      index.set_size(m.cols(), false);

      for (int i=0; i<m.cols(); i++)
	out(i) = min(m.get_col(i), index(i));
    } else {
      out.set_size(m.rows(), false);
      index.set_size(m.rows(), false);

      for (int i=0; i<m.rows(); i++)
	out(i) = min(m.get_row(i), index(i));
    }
    return out;
  }


  //! Return the postion of the maximum element in the vector
  template<class T>
    int max_index(const Vec<T> &in)
  {
    int	maxindex=0;
    for (int i=0;i<in.length();i++)
      if (in[i]>in[maxindex])
	maxindex=i;
    return maxindex;
  }

  //! Return the postion of the maximum element in the matrix
  template<class T>
    void max_index(const Mat<T> &m, int &row, int &col)
  {
    T maxdata = m(0,0);
    int i, j;

    row = col = 0;
    for (i=0; i<m.rows(); i++)
      for (j=0; j<m.cols(); j++)
	if (m(i,j) > maxdata) {
	  row = i;
	  col = j;
	  maxdata = m(i,j);
	}
  }

  //! Return the postion of the minimum element in the vector
  template<class T>
    int min_index(const Vec<T> &in)
  {
    int minindex=0;
    for (int i=0;i<in.length();i++)
      if (in[i]<in[minindex])
	minindex=i;
    return minindex;
  }

  //! Return the postion of the minimum element in the matrix
  template<class T>
    void min_index(const Mat<T> &m, int &row, int &col)
  {
    T mindata = m(0,0);
    int i, j;

    row = col = 0;
    for (i=0; i<m.rows(); i++)
      for (j=0; j<m.cols(); j++)
	if (m(i,j) < mindata) {
	  row = i;
	  col = j;
	  mindata = m(i,j);
	}
  }

  //! The mean value
  double mean(const vec &v);
  //! The mean value
  std::complex<double> mean(const cvec &v);
  //! The mean value
  double mean(const svec &v);
  //! The mean value
  double mean(const ivec &v);
  //! The mean value
  double mean(const mat &m);
  //! The mean value
  std::complex<double> mean(const cmat &m);
  //! The mean value
  double mean(const smat &m);
  //! The mean value
  double mean(const imat &m);

  //! The geometric mean of a vector
  template<class T>
  double geometric_mean(const Vec<T> &v)
  {
    return std::exp(std::log(static_cast<double>(prod(v))) / v.length());
  }

  //! The geometric mean of a matrix
  template<class T>
  double geometric_mean(const Mat<T> &m)
  {
    return std::exp(std::log(static_cast<double>(prod(prod(m)))) 
		    / (m.rows() * m.cols()));
  }

  //! The median
  template<class T>
    double median(const Vec<T> &v)
  {
    Vec<T> invect(v);
    sort(invect);
    return (double)(invect[(invect.length()-1)/2]+invect[invect.length()/2])/2.0;
  }

  //! Calculate the 2-norm: norm(v)=sqrt(sum(abs(v).^2))
  double norm(const cvec &v);

  //! Calculate the 2-norm: norm(v)=sqrt(sum(abs(v).^2))
  template<class T>
  double norm(const Vec<T> &v)
  {
    double E = 0.0;
    for (int i = 0; i < v.size(); i++)
      E += sqr(static_cast<double>(v[i]));

    return std::sqrt(E);
  }

  //! Calculate the p-norm: norm(v,p)=sum(abs(v).^2)^(1/p)
  double norm(const cvec &v, int p);

  //! Calculate the p-norm: norm(v,p)=sum(abs(v).^2)^(1/p)
  template<class T>
  double norm(const Vec<T> &v, int p)
  {
    double E = 0.0;
    for (int i = 0; i < v.size(); i++)
      E += std::pow(fabs(static_cast<double>(v[i])), static_cast<double>(p));

    return std::pow(E, 1.0 / p);
  }

  //! Calculate the frobeniuos norm for s = "fro" (equal to 2-norm)
  double norm(const cvec &v, const std::string &s);

  //! Calculate the frobeniuos norm for s = "fro" (equal to 2-norm)
  template<class T>
  double norm(const Vec<T> &v, const std::string &s)
  {
    it_assert(s == "fro", "norm(): Unrecognised norm");

    double E = 0.0;
    for (int i = 0; i < v.size(); i++)
      E += sqr(static_cast<double>(v[i]));

    return std::sqrt(E);
  }

  /*!
   * Calculate the p-norm of a real matrix
   *
   * p = 1: max(svd(m))
   * p = 2: max(sum(abs(X)))
   *
   * Default if no p is given is the 2-norm
   */
  double norm(const mat &m, int p = 2);

  /*! 
   * Calculate the p-norm of a complex matrix
   *
   * p = 1: max(svd(m))
   * p = 2: max(sum(abs(X)))
   *
   * Default if no p is given is the 2-norm
   */
  double norm(const cmat &m, int p = 2);

  //! Calculate the frobeniuos norm of a matrix for s = "fro"
  double norm(const mat &m, const std::string &s);
  
  //! Calculate the frobeniuos norm of a matrix for s = "fro"
  double norm(const cmat &m, const std::string &s);


  //! The variance of the elements in the vector. Normalized with N-1 to be unbiased.
  double variance(const cvec &v);

  //! The variance of the elements in the vector. Normalized with N-1 to be unbiased.
  template<class T>
    double variance(const Vec<T> &v)
  {
    int len = v.size();
    const T *p=v._data();
    double sum=0.0, sq_sum=0.0;

    for (int i=0; i<len; i++, p++) {
      sum += *p;
      sq_sum += *p * *p;
    }

    return (double)(sq_sum - sum*sum/len) / (len-1);
  }

  //! Calculate the energy: squared 2-norm. energy(v)=sum(abs(v).^2)
  template<class T>
    double energy(const Vec<T> &v)
  {
    return sqr(norm(v));
  }


  //! Return true if the input value \c x is within the tolerance \c tol of the reference value \c xref
  inline bool within_tolerance(double x, double xref, double tol = 1e-14)
  {
    return ( fabs(x-xref) <= tol ) ? true : false; 
  }

  //! Return true if the input value \c x is within the tolerance \c tol of the reference value \c xref
  inline bool within_tolerance(std::complex<double> x, std::complex<double> xref, double tol = 1e-14)
  {
    return ( abs(x-xref) <= tol ) ? true : false; 
  }

  //! Return true if the input vector \c x is elementwise within the tolerance \c tol of the reference vector \c xref
  inline bool within_tolerance(const vec &x, const vec &xref, double tol = 1e-14)
  {
    return ( max(abs(x-xref)) <= tol ) ? true : false; 
  }

  //! Return true if the input vector \c x is elementwise within the tolerance \c tol of the reference vector \c xref
  inline bool within_tolerance(const cvec &x, const cvec &xref, double tol = 1e-14)
  {
    return ( max(abs(x-xref)) <= tol ) ? true : false; 
  }

  //! Return true if the input matrix \c X is elementwise within the tolerance \c tol of the reference matrix \c Xref
  inline bool within_tolerance(const mat &X, const mat &Xref, double tol = 1e-14)
  {
    return ( max(max(abs(X-Xref))) <= tol ) ? true : false; 
  }

  //! Return true if the input matrix \c X is elementwise within the tolerance \c tol of the reference matrix \c Xref
  inline bool within_tolerance(const cmat &X, const cmat &Xref, double tol = 1e-14)
  {
    return ( max(max(abs(X-Xref))) <= tol ) ? true : false; 
  }

  /*!
    \brief Calculate the central moment of vector x
  
    The \f$r\f$th sample central moment of the samples in the vector
    \f$ \mathbf{x} \f$ is defined as
  
    \f[
    m_r = \mathrm{E}[x-\mu]^r = \frac{1}{n} \sum_{i=0}^{n-1} (x_i - \mu)^r
    \f]
    where \f$\mu\f$ is the sample mean.
  */
  double moment(const vec &x, const int r);

  /*!
    \brief Calculate the skewness excess of the input vector x

    The skewness is a measure of the degree of asymmetry of distribution. Negative
    skewness means that the distribution is spread more to the left of the mean than to
    the right, and vice versa if the skewness is positive.

    The skewness of the samples in the vector \f$ \mathbf{x} \f$ is
    \f[
    \gamma_1 = \frac{\mathrm{E}[x-\mu]^3}{\sigma^3}
    \f]
    where \f$\mu\f$ is the mean and \f$\sigma\f$ the standard deviation.

    The skewness is estimated as
    \f[
    \gamma_1 = \frac{k_3}{{k_2}^{3/2}}
    \f]
    where
    \f[
    k_2 = \frac{n}{n-1} m_2
    \f]
    and
    \f[
    k_3 = \frac{n^2}{(n-1)(n-2)} m_3
    \f]
    Here \f$m_2\f$ is the sample variance and \f$m_3\f$ is the 3rd sample
    central moment.
  */
  double skewness(const vec &x);


  /*!
    \brief Calculate the kurtosis excess of the input vector x

    The kurtosis excess is a measure of peakedness of a distribution. 
	The kurtosis excess is defined as
    \f[
    \gamma_2 = \frac{\mathrm{E}[x-\mu]^4}{\sigma^4} - 3
    \f]
    where \f$\mu\f$ is the mean and \f$\sigma\f$ the standard deviation.

    The kurtosis excess is estimated as
    \f[
    \gamma_2 = \frac{k_4}{{k_2}^2}
    \f]
    where
    \f[
    k_2 = \frac{n}{n-1} m_2
    \f]
    and
    \f[
    k_4 = \frac{n^2 [(n+1)m_4 - 3(n-1){m_2}^2]}{(n-1)(n-2)(n-3)}
    \f]
    Here \f$m_2\f$ is the sample variance and \f$m_4\f$ is the 4th sample
    central moment.
  */
  double kurtosisexcess(const vec &x);

    /*!
    \brief Calculate the kurtosis of the input vector x

    The kurtosis is a measure of peakedness of a distribution. The kurtosis 
	is defined as
    \f[
    \gamma_2 = \frac{\mathrm{E}[x-\mu]^4}{\sigma^4}
    \f]
    where \f$\mu\f$ is the mean and \f$\sigma\f$ the standard deviation.
	For a Gaussian variable, the kurtusis is 3.

	See also the definition of kurtosisexcess.
  */
  inline double kurtosis(const vec &x) {return kurtosisexcess(x)+3;}

  //!@}

} // namespace itpp

#endif // #ifndef STAT_H
