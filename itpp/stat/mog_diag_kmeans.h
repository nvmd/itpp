/*!
 * \file 
 * \brief K-means based optimiser for Mixture of Gaussians - header file
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


#ifndef MOG_DIAG_KMEANS_H
#define MOG_DIAG_KMEANS_H

#include <itpp/stat/mog_diag.h>


namespace itpp {
  
  /*! 
    \ingroup MOG
    \brief K-means based optimisation for Mixtures of Gaussians
    \author Conrad Sanderson
    
    This class is an optimiser (trainer) for the parameters of 
    an instance of the MOG_diag class.
    The obtained parameters are typically used as a seed by the 
    MOG_diag_EM class.
  */
  class MOG_diag_kmeans : public MOG_diag {

    public:
    //! Default constructor
    MOG_diag_kmeans() { verbose = false; } 
    
    //! Default destructor
    ~MOG_diag_kmeans() { }
    
    /*!
      \brief Run the k-means algorithm
      
      \param model_in The model to optimise (MOG_diag)
      \param X_in The training data (array of vectors)
      \param max_iter_in Maximum number of iterations. Default is 10.
      \param trust_in The trust factor, where 0 <= \c _trust <= 1.  Default is 0.5.
      \param normalise_in Use normalised distance measure (in effect). Default is true.
      
      \note The higher the trust factor, the more we trust 
      the estimates of covariance matrices and weights.
      Set this to 1.0 only if you have plenty of training data.
      One rule of thumb is 10*D vectors per Gaussian, 
      where D is the dimensionality of the vectors.
      For smaller amounts of data, a lower trust factor
      will help (but not completely avoid) the EM algorithm
      (used in the MOG_diag_em class) from getting stuck
      in a local minimum.
    
      \note Setting \c normalise_in to true causes the the training 
      data to be normalised to zero mean and unit variance prior
      to running the k-means algorithm.  The data is unnormalised
      before returning.  The normalisation helps clustering when 
      the range of values varies greatly between dimensions.
      e.g. dimension 1 may have values in the [-1,+1] interval,
      while dimension 2 may have values in the [-100,+100] interval.
      Without normalisation, the distance between vectors is
      dominated by dimension 2. 
    */
    void run(MOG_diag &model_in, Array<vec> &X_in, int max_iter_in, double trust_in, bool normalise_in);
    
    /*! \brief Enable or disable printing of progress
        \param verbose_in If true, print progress.
    */
    void set_verbose(bool verbose_in) { verbose = verbose_in; }

    protected:
    
    //! squared Euclidean distance between two C vectors
    inline double dist(const double * x, const double * y) const;
    
    void assign_to_means();
    void recalculate_means();
    bool dezombify_means();
    double measure_change() const;
    void initial_means();
    void iterate();
    void calc_means();
    void calc_covs();
    void calc_weights();
    void normalise_vectors();
    void unnormalise_vectors();
    void unnormalise_means();

    //! Maximum number of iterations
    int max_iter;
    
    /*! \brief trust factor, where 0 <= trust <= 1.
        \note The higher the trust factor, the more we trust the estimates of covariance matrices and weights.
     */
    double trust; 
    
    //! Whether we print the progress 
    bool verbose;
    
    //! number of training vectors
    int N;

    //! 'C' pointers to training vectors  
    double ** c_X;
    
    //! means from the previous iteration, used to measure progress
    Array<vec> means_old;
    
    //! 'C' pointers to old means
    double ** c_means_old;
    
    //! contains indices of vectors assigned to each mean
    Array<ivec> partitions;
    
    //! 'C' pointers to partition vectors
    int ** c_partitions;
    
    //! keeps a count of the number of vectors assigned to each mean
    ivec count;
    
    //! 'C' pointer to the count vector
    int * c_count;
    
    private:
    
    vec norm_mu;
    double * c_norm_mu;
    
    vec norm_sd;
    double * c_norm_sd;
    
    vec tmpvec;
    double * c_tmpvec;
      
  
  };

}

//
// functions included in the header file for speed reasons

namespace itpp {

  inline double MOG_diag_kmeans::dist(const double * x, const double * y) const {
    double acc = 0.0;
    for(int d=0;d<D;d++) { double tmp = x[d]-y[d]; acc += tmp*tmp; }
    return(acc);
  }

}

#endif  // #ifndef MOG_DIAG_KMEANS_H

