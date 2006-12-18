/*!
 * \file 
 * \brief generic Mixture of Gaussians (MOG) class - source file
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

#include <itpp/stat/mog_generic.h>
#include <itpp/base/inv.h>
#include <itpp/base/det.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/itfile.h>
#include <itpp/base/itmisc.h>


namespace itpp {


  void MOG_generic::init() { cleanup(); }


  void MOG_generic::init(const int &K_in, const int &D_in, bool full_in) {
    valid = false;
  
    it_assert(K_in >= 0, "MOG_generic::init(): number of Gaussians must be greater than zero");
    it_assert(D_in >= 0, "MOG_generic::init(): dimensionality must be greater than zero");
  
    K = K_in;
    D = D_in;
    full = full_in;
  
    set_means_zero_internal();
    full ? set_full_covs_unity_internal() : set_diag_covs_unity_internal();
    set_weights_uniform_internal();
    setup_misc();
  
    valid = true;
    do_checks = true;
    paranoid = false;

  }


  void MOG_generic::init(Array<vec> &means_in, bool full_in) {
    valid = false;
  
    K = means_in.size();
    D = means_in(0).size();
    full = full_in;
  
    it_assert( check_array_uniformity(means_in), "MOG_generic::init(): 'means' is empty or contains vectors of varying dimensionality" );
    set_means(means_in);
    full ? set_full_covs_unity_internal() : set_diag_covs_unity_internal();
    set_weights_uniform_internal();
    setup_misc();
  
    valid = true;
    do_checks = true;
    paranoid = false;
  }


  void MOG_generic::init(Array<vec> &means_in, Array<vec> &diag_covs_in, vec &weights_in) {
    valid = false;
  
    K = means_in.size();
    D = means_in(0).size();
    full = false;
  
    it_assert( check_array_uniformity(means_in), "MOG_generic::init(): 'means' is empty or contains vectors of varying dimensionality" );
  
    set_means_internal(means_in);
    set_diag_covs_internal(diag_covs_in);
    set_weights_internal(weights_in);
    setup_misc();  
  
    valid = true;
    do_checks = true;
    paranoid = false;
  }


  void MOG_generic::init(Array<vec> &means_in, Array<mat> &full_covs_in, vec &weights_in) {
    valid = false;
  
    K = means_in.size();
    D = means_in(0).size();
    full = true;
  
    it_assert( check_array_uniformity(means_in), "MOG_generic::init(): 'means' is empty or contains vectors of varying dimensionality" );
    set_means_internal(means_in);
    set_full_covs_internal(full_covs_in);
    set_weights_internal(weights_in);
    setup_misc();  
  
    valid = true;
    do_checks = true;
    paranoid = false;
  }

  
  bool MOG_generic::check_size(const vec &x_in) const {
    if(x_in.size() == D) return true;
    return false;
  }
  
  
  bool MOG_generic::check_size(const Array<vec> &X_in) const {
    if(check_array_uniformity(X_in)) return check_size(X_in(0));
    return false;
  }


  bool MOG_generic::check_array_uniformity(const Array<vec> & A) const {
    int rows = A.size();
    int cols = A(0).size();
  
    if(!rows || !cols) return false;
  
    for(int row=1; row<rows; row++)  if (A(row).size() != cols)  return false;
    return true;
  }

  
  void MOG_generic::set_means_zero_internal() {
    means.set_size(K);
    for(int k=0; k<K; k++) { means(k).set_size(D);  means(k) = 0.0; }
    setup_means();
    }


  void MOG_generic::set_means_internal(Array<vec> &means_in) {
    it_assert( (means_in.size() == K), "MOG_generic::set_means_internal(): number of vectors in 'means' is not equivalent to number of Gaussians" );
  
    for(int k=0; k<K; k++)
      it_assert( (means_in(k).size() == D), "MOG_generic::set_means_internal(): dimensionality mismatch between model and one or more vectors in 'means'" );
  
    for(int k=0; k<K; k++)
      for(int d=0; d<D; d++)
	      it_assert( finite(means_in(k)(d)), "MOG_generic::set_means_internal(): 'means' has a non-finite value" );

    means = means_in;
    setup_means();
  }


  void MOG_generic::set_diag_covs_internal(Array<vec> &diag_covs_in) {
    it_assert( (diag_covs_in.size() == K ), "MOG_generic::set_diag_covs_internal(): number of vectors in 'diag_covs' does not match number of Gaussians" );
  
    for(int k=0; k<K; k++)
      it_assert( (diag_covs_in(k).size() == D), "MOG_generic::set_diag_covs_internal(): dimensionality mismatch between model and one or more vectors in 'diag_covs'" );
  
    for(int k=0; k<K; k++)
      for(int d=0; d<D; d++) {
	      it_assert( (diag_covs_in(k)(d) > 0.0), "MOG_generic::set_diag_covs_internal(): 'diag_covs' has a zero or negative value" );
	      it_assert( finite(diag_covs_in(k)(d)), "MOG_generic::set_diag_covs_internal(): 'diag_covs' has a non-finite value" );
      }
      
    full_covs.set_size(0);
    diag_covs = diag_covs_in;
    full = false;
    setup_covs();
  }


  void MOG_generic::set_full_covs_internal(Array<mat> &full_covs_in) {
    it_assert( (full_covs_in.size() == K ), "MOG_generic::set_full_covs_internal(): number of matrices in 'full_covs' does not match number of Gaussians" );
  
    for(int k=0; k<K; k++)
      it_assert( ( (full_covs_in(k).rows() == D) && (full_covs_in(k).cols() == D) ), "MOG_generic::set_full_covs_internal(): dimensionality mismatch between model and one or more matrices in 'full_covs'" );
  
    for(int k=0; k<K; k++)
      for(int i=0; i<D; i++) for(int j=0; j<D; j++) {
	      it_assert( finite(full_covs_in(k)(i,j)), "MOG_generic::set_full_covs_internal(): 'full_covs' has a non-finite value" );
	      if(i==j) it_assert( (full_covs_in(k)(i,j) > 0.0), "MOG_generic::set_full_covs_internal(): 'full_covs' has a zero or negative value on a diagonal" );
	    }

    full_covs = full_covs_in;
    diag_covs.set_size(0);
    full = true;
    setup_covs();
  }


  void MOG_generic::set_weights_internal(vec &weights_in) {
  
    it_assert( (weights_in.size() == K ), "MOG_generic::set_weights_internal(): number of elements in 'weights' does not match number of Gaussians" );
  
    for(int k=0; k<K; k++) {
      it_assert( (weights_in(k) >= 0), "MOG_generic::set_weights_internal(): 'weights' has a negative value" );
      it_assert( finite(weights_in(k)), "MOG_generic::set_weights_internal(): 'weights' has a non-finite value" );
      }
    
    weights = weights_in;
    setup_weights();
  
  }


  void MOG_generic::set_diag_covs_unity_internal() {
    full_covs.set_size(0);
    diag_covs.set_size(K);
  
    for(int k=0; k<K; k++) { diag_covs(k).set_size(D);  diag_covs(k) = 1.0; }
  
    full = false;
    setup_covs();
  }


  void MOG_generic::set_full_covs_unity_internal() {
    full_covs.set_size(K);
    diag_covs.set_size(0);
  
    for(int k=0; k<K; k++) {
      full_covs(k).set_size(D,D);
      full_covs(k) = 0.0;
      for(int d=0;d<D;d++)  full_covs(k)(d,d) = 1.0;
    }
  
    full = true;
    setup_covs();
  }


  void MOG_generic::set_weights_uniform_internal() {
    weights.set_size(K);
    weights = 1.0/K;
    setup_weights();
  }


  void MOG_generic::setup_means() { }
  
  void MOG_generic::setup_covs() {

    double Ddiv2_log_2pi = D/2.0 * std::log(m_2pi);
    log_det_etc.set_size(K);
  
    if(full) {
      full_covs_inv.set_size(K);
      diag_covs_inv_etc.set_size(0);
      for(int k=0;k<K;k++)  full_covs_inv(k) = inv(full_covs(k));
      for(int k=0;k<K;k++)  log_det_etc(k) = -Ddiv2_log_2pi - 0.5*std::log(det(full_covs(k)));
    }
    else {
      full_covs_inv.set_size(0);
      diag_covs_inv_etc.set_size(K);  for(int k=0;k<K;k++) diag_covs_inv_etc(k).set_size(D);    
            
      for(int k=0;k<K;k++) {
        double acc = 0.0;
        vec & diag_cov = diag_covs(k);
        vec & diag_cov_inv_etc = diag_covs_inv_etc(k);
        
        for(int d=0;d<D;d++)  {
          double tmp = diag_cov(d);
          diag_cov_inv_etc(d) = 1.0/(2.0*tmp);
          acc += std::log(tmp);
        }  
        
        log_det_etc(k) = -Ddiv2_log_2pi - 0.5*acc;
      
      }
    }
  }


  void MOG_generic::setup_weights() {
    weights /= sum(weights);
    log_weights = log(weights);
  }


  void MOG_generic::setup_misc() {
    log_max_K = std::log(std::numeric_limits<double>::max() / K);
    tmpvecD.set_size(D);
    tmpvecK.set_size(K);
  }

  
  void MOG_generic::cleanup() {
    
    valid=false;
    do_checks=true;
    K=0;
    D=0;
    
    tmpvecD.set_size(0);
    tmpvecK.set_size(0);
    means.set_size(0);
    diag_covs.set_size(0);
    full_covs.set_size(0);
    weights.set_size(0);
    log_det_etc.set_size(0);
    log_weights.set_size(0);
    full_covs_inv.set_size(0);
    diag_covs_inv_etc.set_size(0);
 
  }
  
  
  void MOG_generic::set_means(Array<vec> &means_in) {
    if(!valid) return;
    set_means_internal(means_in);
  }

  
  void MOG_generic::set_means_zero() {
    if(!valid) return;
    set_means_zero_internal();
  }


  void MOG_generic::set_diag_covs(Array<vec> &diag_covs_in) {
    if(!valid) return;
    set_diag_covs_internal(diag_covs_in);
  }


  void MOG_generic::set_full_covs(Array<mat> &full_covs_in) {
    if(!valid) return;
    set_full_covs_internal(full_covs_in);
  }


  void MOG_generic::set_weights(vec &weights_in) {
    if(!valid) return;
    set_weights_internal(weights_in);
  }


  void MOG_generic::set_diag_covs_unity() {
    if(!valid) return;
    set_diag_covs_unity_internal();
  }


  void MOG_generic::set_full_covs_unity() {
    if(!valid) return;
    set_full_covs_unity_internal();
  }


  void MOG_generic::set_weights_uniform() {
    if(!valid) return;
    set_weights_uniform_internal();
  }


  void MOG_generic::load(const std::string &name_in) {
    valid = false;
  
    it_assert(exist(name_in), "MOG_generic::load(): couldn't access file '"+name_in+"'");
    it_file ff(name_in);
  
    bool contents = ff.exists("means") && ( ff.exists("diag_covs") || ff.exists("full_covs") ) && ff.exists("weights");
    it_assert(contents,"MOG_generic::load(): file '"+name_in+"' doesn't appear to be a model file");
  
    Array<vec> means_in;  ff >> Name("means") >> means_in;
    vec weights_in; ff >> Name("weights") >> weights_in;
  
    if( ff.exists("full_covs") ) {
      Array<mat> full_covs_in; ff >> Name("full_covs") >> full_covs_in;
      init(means_in,full_covs_in,weights_in);
    }
    else {
      Array<vec> diag_covs_in; ff >> Name("diag_covs") >> diag_covs_in;
      init(means_in,diag_covs_in,weights_in);
    }
  
    ff.close();
  
  }


  void MOG_generic::save(const std::string &name_in) const {
    if(!valid) return;
  
    it_file ff(name_in);
  
    ff << Name("means") << means;
    if(full) ff << Name("full_covs") << full_covs;
    else ff << Name("diag_covs") << diag_covs;
    ff << Name("weights") << weights;
  
    ff.close();
  
  }

  void MOG_generic::join(const MOG_generic &B_in) {
  
    if(!valid) return;
    if(!B_in.is_valid()) return;
    
    it_assert( (full == B_in.is_full()), "MOG_generic::join(): given model must be of the same type" );   
    it_assert( (B_in.get_D() == D), "MOG_generic::join(): given model has different dimensionality" );
    it_assert( (B_in.get_K() > 0),  "MOG_generic::join(): given model has no components" );
    
    int new_K = K + B_in.get_K();
    vec new_weights(new_K);
    vec B_in_weights = B_in.get_weights();
        
    double alpha = double(K) / double(new_K);
    double beta = double(B_in.get_K()) / double(new_K);
    
    for(int k=0;k<K;k++)  new_weights(k) = alpha * weights(k);
    for(int k=K;k<new_K;k++)  new_weights(k) = beta * B_in_weights(k);
        
    Array<vec> new_means = concat( means, B_in.get_means() );
    
    if(full) {
      Array<mat> new_full_covs = concat(full_covs, B_in.get_full_covs());
      init(new_means, new_full_covs, new_weights);
    }
    else {
      Array<vec> new_diag_covs = concat(diag_covs, B_in.get_diag_covs());
      init(new_means, new_diag_covs, new_weights);
    }
  }


  void MOG_generic::convert_to_diag_internal() {
    if(!full) return;
    
    diag_covs.set_size(K);
    for(int k=0;k<K;k++)  diag_covs(k) = diag(full_covs(k));
    full_covs.set_size(0);
  
    full = false;
    setup_covs();
  }


  void MOG_generic::convert_to_diag() {
    if(!valid) return;
    convert_to_diag_internal();
  }
  
  
  void MOG_generic::convert_to_full_internal() {
    if(full) return;
  
    full_covs.set_size(K);
    for(int k=0;k<K;k++)  full_covs(k) = diag(diag_covs(k));
    diag_covs.set_size(0);
  
    full = true;
    setup_covs();
  }

  void MOG_generic::convert_to_full() {
    if(!valid) return;
    convert_to_full_internal();
  }
  
} // namespace itpp
