/*!
 * \file 
 * \brief Generic Mixture of Gaussians (MOG) class - header file
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

#ifndef MOG_GENERIC_H
#define MOG_GENERIC_H

#include <itpp/base/mat.h>
#include <itpp/base/array.h>
#include <itpp/base/logexpfunc.h>


namespace itpp {

  /*! 
    \ingroup MOG
    \brief Generic Mixture of Gaussians (MOG) class. Used as a base 
		for other MOG classes. 
    \author Conrad Sanderson
    
    Used for representing a statistical distribution as a 
    convex combination of multi-variate Gaussian functions.
    Also known as a Gaussian Mixture Model.
    This class handles both full and diagonal covariance matrices,
    allows loading and saving of the MOG's parameters,
    as well as calculation of likelihoods. The parameters 
    are set by the user or an optimisation algorithm
    (for example, see the MOG_diag_EM class).
    For speed and space reasons, diagonal and full
    covariance matrices are stored and handled separately. 
  */
  class MOG_generic {
  
  public:
    
    /*! \brief Default constructor
        \note An empty model is created.
              The likelihood functions are not useable
              until the model's parameters are set 
    */
    MOG_generic() { init(); }
    
    /*! \brief Construct the MOG_generic object by loading the parameters from a model file
        \param _name The model's filename
    */ 
    MOG_generic(const std::string &_name) { load(_name); }

    /*! \brief construct a default model (all Gaussians have zero mean and unit variance for all dimensions)
        \param _K Number of Gaussians
        \param _D Dimensionality
        \param _full If true, use full covariance matrices; if false, use diagonal covariance matrices. Default = false. 
    */
    MOG_generic(int _K, int _D, bool _full=false) { init(_K, _D, _full); }

    /*! \brief Construct a model using user supplied mean vectors
        \param _means Array of mean vectors
        \param _full If true, use full covariance matrices; if false, use diagonal covariance matrices. Default = false.
        \note  The number of mean vectors specifies the number of Gaussians.
        The covariance matrices are set to the identity matrix.
        The weights for all Gaussians are the same, equal to 1/K, where K is the number of Gaussians
    */
    MOG_generic(Array<vec> &_means, bool _full=false) { init(_means, _full); }

    /*! \brief Construct a model using user supplied parameters (diagonal covariance version) 
        \param _means Array of mean vectors
        \param _diag_covs Array of vectors representing diagonal covariances
        \param _weights vector of weights
        \note  The number of mean vectors, covariance vectors and weights must be the same
    */
    MOG_generic(Array<vec> &_means, Array<vec> &_diag_covs, vec &_weights) { init(_means, _diag_covs, _weights); }

    /*! \brief Construct a model using user supplied parameters (full covariance version) 
        \param _means Array of mean vectors
        \param _diag_covs Array of full covariance matrices
        \param _weights vector of weights
        \note  The number of mean vectors, covariance matrices and weights must be the same
    */
    MOG_generic(Array<vec> &_means, Array<mat> &_full_covs, vec &_weights) { init(_means, _full_covs, _weights); }

    //! Default destructor
    virtual ~MOG_generic() { cleanup(); }
  
    /*! \brief Initialise the model to be empty.
        \note The likelihood functions are not useable
              until the model's parameters are set 
    */
    void init();

    /*! \brief initialise the model so that all Gaussians have zero mean and unit variance for all dimensions  
        \param _K Number of Gaussians
        \param _D Dimensionality
        \param _full If true, use full covariance matrices; if false, use diagonal covariance matrices. Default = false. 
    */
    void init(int _K, int _D, bool _full);

    /*! \brief Initialise the model using user supplied mean vectors
        \param _means Array of mean vectors
        \param _full If true, use full covariance matrices; if false, use diagonal covariance matrices. Default = false.
        \note  The number of mean vectors specifies the number of Gaussians.
        The covariance matrices are set to the identity matrix.
        The weights for all Gaussians are the same, equal to 1/K, where K is the number of Gaussians
    */
    void init(Array<vec> &_means, bool _full);

    /*! \brief Initialise the model using user supplied parameters (diagonal covariance version) 
        \param _means Array of mean vectors
        \param _diag_covs Array of vectors representing diagonal covariances
        \param _weights vector of weights
        \note  The number of mean vectors, covariance vectors and weights must be the same
      */
    void init(Array<vec> &_means, Array<vec> &_diag_covs, vec &_weights);

    /*! \brief Initialise the model using user supplied parameters (full covariance version) 
        \param _means Array of mean vectors
        \param _diag_covs Array of covariance matrices
        \param _weights vector of weights
        \note  The number of mean vectors, covariance matrices and weights must be the same
    */
    void init(Array<vec> &_means, Array<mat> &_full_covs, vec &_weights);

    /*! \brief Release memory used by the model. The model will be empty.
        \note The likelihood functions are not useable
              until the model's parameters are re-initialised
    */
    virtual void cleanup();
  
    //! Returns true if the model's parameters are valid
    bool is_valid() const { return valid; }

    //! Returns true if the model has full covariance matrices
    bool is_full() const { return full; }

    //! Return the number of Gaussians 
    int get_K() const { if(valid) return(K); else return(0); }
    
    //! Return the dimensionality
    int get_D() const { if(valid) return(D); else return(0); }
  
    //! Obtain a copy of the weight vector
    vec get_weights() const { vec tmp;  if(valid) { tmp = weights; } return tmp; }

    //! Obtain a copy of the array of mean vectors
    Array<vec> get_means() const { Array<vec> tmp; if(valid) { tmp = means; } return tmp; }

    //! Obtain a copy of the array of diagonal covariance vectors
    Array<vec> get_diag_covs() const { Array<vec> tmp; if(valid && !full) { tmp = diag_covs; } return tmp; }

    //! Obtain a copy of the array of full covariance matrices 
    Array<mat> get_full_covs() const { Array<mat> tmp; if(valid && full) { tmp = full_covs; } return tmp; }
  
    /*! \brief Set the means of the model
        \note The number of means must match the number of Gaussians in the model
    */ 
    void set_means(Array<vec> &_means);

    /*! \brief Set the diagonal covariance vectors of the model
        \note  The number of diagonal covariance vectors must match the number of Gaussians in the model
    */
    void set_diag_covs(Array<vec> &_diag_covs);

    /*! \brief Set the full covariance matrices of the model
        \note  The number of covariance matrices must match the number of Gaussians in the model
    */
    void set_full_covs(Array<mat> &_full_covs);

    /*! \brief Set the weight vector of the model
        \note  The number of elements in the weight vector must match the number of Gaussians in the model
    */
    void set_weights(vec &_weights);
  
    //! Set the means in the model to be zero
    void set_means_zero();

    //! Set the diagonal covariance vectors to be unity
    void set_diag_covs_unity();

    //! Set the full covariance matrices to be unity
    void set_full_covs_unity();

    //! Set all the weights to 1/K, where K is the number of Gaussians 
    void set_weights_uniform();
  
    /*! \brief Enable/disable internal checks for likelihood functions
        \param _do_checks If true, checks are enabled; if false, checks are disabled
        \note Disabling checks will provide a speedup in the likelihood functions.
              Disable them only when you're happy that everything is working correctly.
    */
    void set_checks(bool _do_checks) { do_checks = _do_checks; }
  
    /*! \brief Enable/disable paranoia about numerical stability
        \param _paranoid If true, calculate likelihoods using a safer, but slower method.
    */
    void set_paranoid(bool _paranoid) { paranoid = _paranoid; }
    
    /*! \brief Initialise the model by loading the parameters from a model file
        \param _name The model's filename
    */
    virtual void load(const std::string &_name);

    /*! \brief Save the model's parameters to a model file
        \param _name The model's filename
    */
    virtual void save(const std::string &_name) const;

    /*! \brief Mathematically join the model with a user supplied model
        \param _B user supplied model
        \note The Arrays of mean vectors and covariance vectors/matrices from the two models
              are simply concatenated, while the weights of the resultant model are a function 
              of the original weights and numbers of Gaussians from both models.
              Specifically,
              \f$ w_{new} = [ \alpha \cdot w_{A} ~~~ \beta \cdot w_{B} ]^T \f$,
              where \f$ w_{new} \f$ is the new weight vector,
              \f$ w_{A} \f$ and \f$ w_{B} \f$ are the weight vectors from model A and B,
              while \f$ \alpha = K_A / (K_A + K_B) \f$
              and \f$ \beta = 1-\alpha \f$.
              In turn, \f$ K_A \f$ and \f$ K_B \f$ are the numbers of Gaussians in model A and B, respectively.   
                
              See <a href="http://dx.doi.org/10.1016/j.patcog.2005.07.001">On transforming statistical models...</a>
              for more information.
    */
    virtual void join(const MOG_generic &_B);
  
    /*! \brief Convert the model to use diagonal covariances
                
        \note If the model is already diagonal, nothing is done. 
              If the model has full covariance matrices,
              this results in irreversible information loss
              (in effect the off-diagonal covariance elements are now zero)
    */ 
    virtual void convert_to_diag();
    
    /*! \brief Convert the model to have full covariance matrices
        \note If the model has full covariance matrices, nothing is done.
              If the model has diagonal covariances, the off-diagonal 
              elements in the full covariance matrices are set to zero.
    */
    virtual void convert_to_full();
  
    //! calculate the log likelihood of vector \c _x using only Gaussian \c k 
    virtual inline double log_lhood_single_gaus(const vec &_x, const int k);
    
    //! calculate the log likelihood of vector \c _x 
    virtual inline double log_lhood(const vec &_x);

    //! calculate the likelihood of vector \c _x 
    virtual inline double lhood(const vec &_x);
    
    //! calculate the average log likelihood of an array of vectors \c _X
    virtual inline double avg_log_lhood(const Array<vec> &_X);

  protected:
    
    //! indicates whether checks on input data are done
    bool do_checks;
    
    //! indicates whether the parameters are valid
    bool valid;
    
    //! indicates whether we are using full or diagonal covariance matrices
    bool full;
  
    //! indicates whether we are paranoid about numerical stability 
    bool paranoid;
    
    //! number of gaussians
    int K;
    
    //! dimensionality
    int D;
    
    //! means
    Array<vec> means;
    
    //! diagonal covariance matrices, stored as vectors
    Array<vec> diag_covs;
    
    //! full covariance matrices
    Array<mat> full_covs;
    
    //! weights
    vec weights;
  
    //! Pre-calcualted std::log(std::numeric_limits<double>::max() / K), where K is the number of Gaussians
    double log_max_K;

    /*! \brief Gaussian specific pre-calcualted constants
        \note  Vector of pre-calculated \f$ -\frac{D}{2}\log(2\pi) -\frac{1}{2}\log(|\Sigma|) \f$ for each Gaussian,
               where \f$ D \f$ is the dimensionality and \f$ |\Sigma| \f$ is the determinant for the Gaussian's
               covariance matrix \f$ \Sigma \f$.
    */    
    vec log_det_etc;

    //! Pre-calculated log versions of the weights
    vec log_weights;

    //! Pre-calcuated inverted version of each full covariance matrix 
    Array<mat> full_covs_inv;

    //! Pre-calcuated inverted version of each diagonal covariance vector, where the covariance elements are first multiplied by two 
    Array<vec> diag_covs_inv_etc;
    
    //! Check if vector \c _x has the same dimensionality as the model
    bool check_size(const vec &_x) const;

    //! Check if all vectors in Array \c _X have the same dimensionality as the model
    bool check_size(const Array<vec> &_X) const;
    
    //! Check if all vectors in Array \c _X have the same dimensionality
    bool check_array_uniformity(const Array<vec> & A) const;

    void set_means_internal(Array<vec> &_means);
    void set_diag_covs_internal(Array<vec> &_diag_covs);
    void set_full_covs_internal(Array<mat> &_full_covs);
    void set_weights_internal(vec &_weigths);
  
    void set_means_zero_internal();
    void set_diag_covs_unity_internal();
    void set_full_covs_unity_internal();
    void set_weights_uniform_internal();

    void convert_to_diag_internal();
    void convert_to_full_internal();

    //! additional processing of mean vectors, done as the last step of mean initialisation
    virtual void setup_means();

    //! additional processing of covariance vectors/matrices, done as the last step of covariance initialisation
    virtual void setup_covs();

    //! additional processing of the weight vector, done as the last step of weight initialisation
    virtual void setup_weights();

    //! additional processing of miscellaneous parameters, done as the last step of overall initialisation
    virtual void setup_misc();
  
    virtual inline double log_lhood_single_gaus_internal(const vec &_x, const int k);
    virtual inline double log_lhood_internal(const vec &_x);
    virtual inline double lhood_internal(const vec &_x);
  
  private:
    vec tmpvecD;
    vec tmpvecK;

  };
  
} // namespace itpp


//
// functions included in the header file for speed reasons

namespace itpp {

  inline double MOG_generic::log_lhood_single_gaus_internal(const vec &_x, const int k) {
   
    const vec & mean = means(k);
   
    if(full) {
      for(int d=0;d<D;d++)  tmpvecD[d] = _x[d] - mean[d];
      double tmpval = tmpvecD*(full_covs_inv(k)*tmpvecD);
      return(log_det_etc[k] - 0.5*tmpval);
    }
    else {
      const vec & diag_cov_inv_etc = diag_covs_inv_etc(k);
    
      double acc = 0.0;
    
      for(int d=0; d<D; d++) {
        double tmpval = _x[d] - mean[d];
        acc += (tmpval*tmpval) * diag_cov_inv_etc[d];
      }
      return(log_det_etc[k] - acc);
    }
  
  }
  
  inline double MOG_generic::log_lhood_single_gaus(const vec &_x, const int k) {
    if(do_checks) {
      it_assert(valid, "MOG_generic::log_lhood_single_gaus(): model not valid");
      it_assert(check_size(_x), "MOG_generic::log_lhood_single_gaus(): x has wrong dimensionality");
      it_assert( ( (k>=0) && (k<K) ), "MOG_generic::log_lhood_single_gaus(): k specifies a non-existant Gaussian");
    }
    return log_lhood_single_gaus_internal(_x,k);
  }

  
  inline double MOG_generic::log_lhood_internal(const vec &_x) {
    
    bool danger = paranoid;
    
    for(int k=0;k<K;k++)  {
      double tmp = log_weights[k] + log_lhood_single_gaus_internal(_x,k); 
      tmpvecK[k] = tmp;
      
      if(tmp >= log_max_K)  danger = true;
    }
    
    if(danger) {
      double log_sum = tmpvecK[0];  for(int k=1; k<K; k++)  log_sum = log_add( log_sum, tmpvecK[k] );
      return(log_sum);
    }
    else {
      double sum = 0.0; for(int k=0;k<K;k++) sum += std::exp(tmpvecK[k]);
      return(std::log(sum));
    }
  }
  
  
  inline double MOG_generic::log_lhood(const vec &_x) {
    if(do_checks) {
      it_assert(valid, "MOG_generic::log_lhood(): model not valid");
      it_assert(check_size(_x), "MOG_generic::log_lhood(): x has wrong dimensionality");
    }
    return log_lhood_internal(_x);
  }

        
  inline double MOG_generic::lhood_internal(const vec &_x) {
    
    bool danger = paranoid;
    
    for(int k=0;k<K;k++)  {
      double tmp = log_weights[k] + log_lhood_single_gaus_internal(_x,k); 
      tmpvecK[k] = tmp;
      
      if(tmp >= log_max_K)  danger = true;
    }
    
    if(danger) {
      double log_sum = tmpvecK[0];  for(int k=1; k<K; k++)  log_sum = log_add( log_sum, tmpvecK[k] );
      return(trunc_exp(log_sum));
    }
    else {
      double sum = 0.0; for(int k=0;k<K;k++) sum += std::exp(tmpvecK[k]);
      return(sum);
    }
  }
        
        
  inline double MOG_generic::lhood(const vec &_x) {
    
    if(do_checks) {
      it_assert(valid, "MOG_generic::lhood(): model not valid");
      it_assert(check_size(_x), "MOG_generic::lhood(): x has wrong dimensionality");
    }
    return lhood_internal(_x);
  }
  

  inline double MOG_generic::avg_log_lhood(const Array<vec> &_X) {
    if(do_checks) {
      it_assert(valid, "MOG_generic::avg_log_lhood(): model not valid");
      it_assert(check_size(_X), "MOG_generic::avg_log_lhood(): X is empty or at least one vector has the wrong dimensionality");
    }
  
    const int N = _X.size();
    double acc = 0.0;
    for(int n=0;n<N;n++)  acc += log_lhood_internal(_X(n));
    return(acc/N);
  }

} // namespace itpp

#endif // #ifndef MOG_GENERIC_H
