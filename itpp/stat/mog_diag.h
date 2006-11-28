/*!
 * \file 
 * \brief Diagonal Mixture of Gaussians class - header file
 * \author Conrad Sanderson
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

#ifndef MOG_DIAG_H
#define MOG_DIAG_H

#include <itpp/stat/mog_generic.h>


namespace itpp {

  /*! 
    \ingroup MOG
    \brief Diagonal Mixture of Gaussians (MOG) class 
    \author Conrad Sanderson
    
    Used for representing a statistical distribution as a 
    convex combination of multi-variate Gaussian functions.
    Also known as a Gaussian Mixture Model.
    This class allows loading and saving of the MOG's parameters,
    as well as calculation of likelihoods. The parameters 
    are set by the user or an optimisation algorithm
    (for example, see the MOG_diag_EM class).
    
    \note This class is optimised for diagonal covariance matrices.
          For speed reasons it uses C style arrays for direct access to memory. 
  */
  class MOG_diag : public MOG_generic {

    public:
    
    /*! \brief Default constructor
        \note An empty model is created.
              The likelihood functions are not useable
              until the model's parameters are set 
    */
    MOG_diag() { zero_all_ptrs(); init(); }

    /*! \brief Construct the MOG_diag object by loading the parameters from a model file
        \param _name The model's filename
    */ 
    MOG_diag(const std::string &name) { zero_all_ptrs(); load(name); }

    /*! \brief construct a default model (all Gaussians have zero mean and unit variance for all dimensions)
        \param _K Number of Gaussians
        \param _D Dimensionality
        \param _full Ignored.  Present for compatability with the MOG_generic class 
    */
    MOG_diag(int _K, int _D, bool _full=false) { zero_all_ptrs(); init(_K,_D,false); }

    /*! \brief Construct a model using user supplied mean vectors
        \param _means Array of mean vectors
        \param _full Ignored.  Present for compatability with the MOG_generic class
        \note  The number of mean vectors specifies the number of Gaussians.
        The covariance matrices are in effect set equal to the identity matrix.
        The weights for all Gaussians are the same, equal to 1/K, where K is the number of Gaussians
    */
    MOG_diag(Array<vec> &_means, bool _full=false) { zero_all_ptrs(); init(_means,false);  }

    /*! \brief Construct a model using user supplied parameters (diagonal covariance version) 
        \param _means Array of mean vectors
        \param _diag_covs Array of vectors representing diagonal covariances
        \param _weights vector of weights
        \note  The number of mean vectors, covariance vectors and weights must be the same
    */
    MOG_diag(Array<vec> &_means, Array<vec> &_diag_covs, vec &_weights) { zero_all_ptrs(); init(_means,_diag_covs,_weights); }

    /*! \brief Construct a model using user supplied parameters (full covariance version) 
        \param _means Array of mean vectors
        \param _diag_covs Array of full covariance matrices
        \param _weights vector of weights
        \note  The full covariance matrices are converted to be diagonal.
               The number of mean vectors, covariance matrices and weights must be the same.
    */
    MOG_diag(Array<vec> &_means, Array<mat> &_full_covs, vec &_weights) { zero_all_ptrs(); init(_means,_full_covs,_weights); convert_to_diag(); }

    //! Default destructor
    ~MOG_diag() { cleanup(); }
    
    /*! \brief Release memory used by the model. The model will be empty.
        \note The likelihood functions are not useable
              until the model's parameters are re-initialised
    */
    void cleanup() { free_all_ptrs(); MOG_generic::cleanup(); }

    /*! \brief Initialise the model by loading the parameters from a model file.
        \param _name The model's filename
        \note If the model file contains a full covariance matrix model,
              the covariance matrices will be converted to be diagonal after loading. 
    */
    void load(const std::string &_name);
    
    //! Do nothing.  Present for compatability with the MOG_generic class.
    void convert_to_full() {};

    //! calculate the log likelihood of C vector \c c_x using only Gaussian \c k 
    inline double log_lhood_single_gaus(const double * c_x, const int k) const;

    //! calculate the log likelihood of IT++ vector \c _x using only Gaussian \c k 
    inline double log_lhood_single_gaus(const vec &_x, const int k) const;

    //! calculate the log likelihood of C vector \c c_x 
    inline double log_lhood(const double * c_x);

    //! calculate the log likelihood of IT++ vector \c _x 
    inline double log_lhood(const vec &_x);

    //! calculate the likelihood of C vector \c c_x 
    inline double lhood(const double * c_x);

    //! calculate the likelihood of IT++ vector \c _x 
    inline double lhood(const vec &_x);

    //! calculate the average log likelihood of an array of C vectors ( \c c_X ) 
    inline double avg_log_lhood(const double ** c_X, int N);

    //! calculate the average log likelihood of an array of IT++ vectors ( \c _X ) 
    inline double avg_log_lhood(const Array<vec> & _X); 
    
    protected:
    
    void setup_means();
    void setup_covs();
    void setup_weights();
    void setup_misc();
    
    inline double log_lhood_single_gaus_internal(const double * c_x, const int k) const;
    inline double log_lhood_single_gaus_internal(const vec &_x, const int k) const;
    inline double log_lhood_internal(const double * c_x);
    inline double log_lhood_internal(const vec &_x);
    inline double lhood_internal(const double * c_x);
    inline double lhood_internal(const vec &_x);

    //! Enable C style access to an Array of vectors (vec)
    double ** enable_c_access(Array<vec> & _A);

    //! Enable C style access to an Array of vectors (ivec)
    int ** enable_c_access(Array<ivec> & _A);

    //! Enable C style access to a vector (vec)
    double * enable_c_access(vec & _v);

    //! Enable C style access to a vector (ivec)
    int * enable_c_access(ivec & _v);
    
    //! Disable C style access to an Array of vectors (vec)
    double ** disable_c_access(double ** A);

    //! Disable C style access to an Array of vectors (ivec)
    int ** disable_c_access(int ** A);

    //! Disable C style access to a vector (vec)
    double * disable_c_access(double * v);

    //! Disable C style access to a vector (ivec)
    int * disable_c_access(int * v);
    
    void zero_all_ptrs();
    void free_all_ptrs();
    
    //! pointers to the mean vectors
    double ** c_means;

    //! pointers to the covariance vectors
    double ** c_diag_covs;
    
    //! pointers to the inverted covariance vectors
    double ** c_diag_covs_inv_etc;

    //! pointer to the weight vector
    double * c_weights;  
    
    //! pointer to the log version of the weight vector
    double * c_log_weights;

    //! pointer to the log_det_etc vector
    double * c_log_det_etc;

    private:
    
    vec tmpvecK;
    double * c_tmpvecK;

  };

}
  
//
// functions included in the header file for speed reasons

namespace itpp {

  inline double MOG_diag::log_lhood_single_gaus_internal(const double * c_x, const int k) const {
  
    const double * c_mean = c_means[k];
    const double * c_diag_cov_inv_etc = c_diag_covs_inv_etc[k];  
  
    double acc = 0.0;
  
    for(int d=0; d<D; d++) {
      double tmp_val = c_x[d] - c_mean[d];
      acc += (tmp_val*tmp_val) * c_diag_cov_inv_etc[d];
    }
    return(c_log_det_etc[k] - acc);
  }


  inline double MOG_diag::log_lhood_single_gaus_internal(const vec &_x, const int k) const {
    return log_lhood_single_gaus_internal(_x._data(), k);
  }


  inline double MOG_diag::log_lhood_single_gaus(const double * c_x, const int k) const {
    if(do_checks) {
      it_assert(valid, "MOG_diag::log_lhood_single_gaus(): model not valid");
      it_assert( ( (k>=0) && (k<K) ), "MOG::log_lhood_single_gaus(): k specifies a non-existant Gaussian");
    }
    return log_lhood_single_gaus_internal(c_x,k);
  }


  inline double MOG_diag::log_lhood_single_gaus(const vec &_x, const int k) const {
    if(do_checks) {
      it_assert(valid, "MOG_diag::log_lhood_single_gaus(): model not valid");
      it_assert(check_size(_x), "MOG_diag::log_lhood_single_gaus(): x has wrong dimensionality");
      it_assert( ( (k>=0) && (k<K) ), "MOG::log_lhood_single_gaus(): k specifies a non-existant Gaussian");
    }
    return log_lhood_single_gaus_internal(_x._data(),k);
  }


  inline double MOG_diag::log_lhood_internal(const double * c_x) {
    
    bool danger = paranoid;

    for(int k=0;k<K;k++)  {
      double tmp = c_log_weights[k] + log_lhood_single_gaus_internal(c_x,k); 
      c_tmpvecK[k] = tmp;
      
      if(tmp >= log_max_K)  danger = true;
    }
  
    
    if(danger) {
      double log_sum = c_tmpvecK[0];  for(int k=1; k<K; k++)  log_sum = log_add( log_sum, c_tmpvecK[k] );
      return(log_sum);
    }
    else {
      double sum = 0.0; for(int k=0;k<K;k++) sum += std::exp(c_tmpvecK[k]);
      return(std::log(sum));
    }
  }
  
  
  inline double MOG_diag::log_lhood_internal(const vec &_x)  {
    return log_lhood_internal(_x._data());
  }


  inline double MOG_diag::log_lhood(const vec &_x) {
    if(do_checks) {
      it_assert(valid, "MOG_diag::log_lhood(): model not valid");
      it_assert(check_size(_x), "MOG_diag::log_lhood(): x has wrong dimensionality");
    }
    return log_lhood_internal(_x._data());
  }


  inline double MOG_diag::log_lhood(const double * c_x) {
    if(do_checks) {
      it_assert(valid, "MOG_diag::log_lhood(): model not valid");
      it_assert( (c_x != 0), "MOG_diag::log_lhood(): c_x is a null pointer");
    }
  
    return log_lhood_internal(c_x);
  }

        
  inline double MOG_diag::lhood_internal(const double * c_x) {
    
    bool danger = paranoid;

    for(int k=0;k<K;k++)  {
      double tmp = c_log_weights[k] + log_lhood_single_gaus_internal(c_x,k); 
      c_tmpvecK[k] = tmp;
      
      if(tmp >= log_max_K)  danger = true;
    }
  
    
    if(danger) {
      double log_sum = c_tmpvecK[0];  for(int k=1; k<K; k++)  log_sum = log_add( log_sum, c_tmpvecK[k] );
      return(trunc_exp(log_sum));
    }
    else {
      double sum = 0.0; for(int k=0;k<K;k++) sum += std::exp(c_tmpvecK[k]);
      return(sum);
    }
  }
  
  inline double MOG_diag::lhood_internal(const vec &_x) { return lhood_internal(_x._data()); }
  
  inline double MOG_diag::lhood(const vec &_x) {
    if(do_checks) {
      it_assert(valid, "MOG_diag::lhood(): model not valid");
      it_assert(check_size(_x), "MOG_diag::lhood(): x has wrong dimensionality");
    }
    return lhood_internal(_x._data());
  }


  inline double MOG_diag::lhood(const double * c_x) {
    if(do_checks) {
      it_assert(valid, "MOG_diag::lhood(): model not valid");
      it_assert( (c_x != 0), "MOG_diag::lhood(): c_x is a null pointer");
    }
  
    return lhood_internal(c_x);
  }


  inline double MOG_diag::avg_log_lhood(const double ** c_X, const int N) {
    if(do_checks) {
      it_assert(valid, "MOG_diag::avg_log_lhood(): model not valid");
      it_assert( (c_X != 0), "MOG_diag::avg_log_lhood(): c_X is a null pointer");
      it_assert( (N >= 0), "MOG_diag::avg_log_lhood(): N is zero or negative");
    }
      
    double acc = 0.0;  for(int n=0;n<N;n++) acc += log_lhood_internal(c_X[n]);
    return(acc/N);
  }


  inline double MOG_diag::avg_log_lhood(const Array<vec> &_X) {
    if(do_checks) {
      it_assert(valid, "MOG_diag::avg_log_lhood(): model not valid");
      it_assert(check_size(_X), "MOG_diag::avg_log_lhood(): X is empty or at least one vector has the wrong dimensionality");
    } 
    const int N = _X.size();
    double acc = 0.0;
    for(int n=0;n<N;n++)  acc += log_lhood_internal(_X(n)._data());
    return(acc/N);
  }
  
}
  
#endif // #ifndef MOG_DIAG_H

