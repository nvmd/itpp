/*!
 * \file 
 * \brief Expectation Maximisation (EM) based optimisers for MOG - header file
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

#ifndef MOG_DIAG_EM_H
#define MOG_DIAG_EM_H

#include <itpp/stat/mog_diag.h>


namespace itpp {

  /*! 
    \ingroup MOG
    \brief Expectation Maximisation (EM) based optimisers for Mixtures 
		of Gaussians
    \author Conrad Sanderson
    
		This class is an optimiser (trainer) for the parameters of an instance of
		the MOG_diag class. The seed values (starting points) are typically
		obtained via the MOG_diag_kmeans class.  See 
		<a href="http://ieeexplore.ieee.org/xpl/abs_free.jsp?arNumber=1561601">
		[CSB06]</a> and the references therein for detailed mathematical
		descriptions.

		- [CSB06] F. Cardinaux, C. Sanderson and S. Bengio, "User authentication
		via adapted statistical models of face images", IEEE Transactions on
		Signal Processing, Vol 54, pp. 361- 373, January 2006
  */
  class MOG_diag_EM : public MOG_diag {

    public:

    //! Default constructor
    MOG_diag_EM() { verbose=false; }

    //! Default destructor
    ~MOG_diag_EM() { }
  
   
    /*!
      \brief Run the Maximum Likelihood (ML) version of the EM algorithm
      
      \param _model The model to optimise (MOG_diag)
      \param _X The training data (array of vectors)
      \param _max_iter Maximum number of iterations. Default is 10.
      \param _var_floor Variance floor (lowest allowable variance).  Default is 0.0 (but see the note below)
      \param _weight_floor  Weight floor (lowest allowable weight).  Default is 0.0 (but see the note below)
      
      \note The variance and weight floors are set to std::numeric_limits<double>::min()
            if they are below that value. As such, they are machine dependant.
            The largest allowable weight floor is 1/K, where K is the number of Gaussians.       
    */
    void ml(MOG_diag &_model, Array<vec> &_X, int _max_iter, double _var_floor, double _weight_floor);

    /*!
      \brief NOT YET IMPLEMENTED. Run the Maximum a Posteriori (MAP) version of the EM algorithm.
      
      \param _model The model to optimise (MOG_diag)
      \param _model The model representing the prior
      \param _X The training data (array of vectors)
      \param _max_iter Maximum number of iterations
      \param _alpha Coefficient for combining the parameters with the prior.  0 <= _alpha <= 1.
      \param _var_floor Variance floor (lowest allowable variance).  Set to 0.0 to use the default.
      \param _weight_floor  Weight floor (lowest allowable weight).  Set to 0.0 to use the default.
      
      \note NOT YET IMPLEMENTED. 
      \note The variance and weight floors are set to std::numeric_limits<double>::min()
            if they are below that value.
            The largest allowable weight floor is 1/K, where K is the number of Gaussians.       
    */
    void map(MOG_diag &_out_model, MOG_diag &_prior_model, Array<vec> &_X, int _max_iter, double _alpha, double _var_floor, double _weight_floor);

    /*! \brief Enable or disable printing of progress
        \param _verbose If true, print progress.
    */
    void set_verbose(bool _verbose) { verbose = _verbose; }

    protected:
  
    //! Whether we print the progress 
    bool verbose;
  
    //! number of training vectors
    int N;

    //! Maximum number of iterations
    int max_iter;
  
    //! 'C' pointers to training vectors
    double ** c_X;
  
    double var_floor;
    double weight_floor;
  
    void inline update_internals();
    void inline sanitise_params();
    double ml_update_params();
    void ml_iterate();
  
    private:
    
    vec tmpvecK;
    vec tmpvecD;  
    vec acc_loglhood_K;  
    
    Array<vec> acc_means;
    Array<vec> acc_covs;
    
    double * c_tmpvecK;
    double * c_tmpvecD;
    double * c_acc_loglhood_K;
    
    double ** c_acc_means;
    double ** c_acc_covs;
  
    
  };

}

#endif // #ifndef MOG_DIAG_EM_H

