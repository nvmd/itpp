/*!
 * \file 
 * \brief kmeans based optimiser for Mixture of Gaussians - source file
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


#include <itpp/stat/mog_diag_kmeans.h>

#include <iostream>


namespace itpp {


  void MOG_diag_kmeans::assign_to_means() {

    for(int k=0;k<K;k++) c_count[k] = 0;
  
    for(int n=0;n<N;n++) {
    
      int k = 0;
      double min_dist = dist( c_means[k], c_X[n] );
      int k_winner = k;
    
      for(int k=1;k<K;k++) {
        double tmp_dist = dist( c_means[k], c_X[n] );
        if(tmp_dist < min_dist) { min_dist = tmp_dist; k_winner = k; }
      }
    
      c_partitions[ k_winner ][ count[k_winner] ] = n;
      c_count[k_winner]++;
    }
  }


  void MOG_diag_kmeans::recalculate_means() {

    for(int k=0;k<K;k++) {

      for(int d=0;d<D;d++)  c_tmpvec[d] = 0.0;

      int Nk = c_count[k];
    
      for(int n=0;n<Nk;n++) {
        double * x = c_X[ c_partitions[k][n] ];
        for(int d=0;d<D;d++)  c_tmpvec[d] += x[d];
      }
    
      if(Nk > 0) {
        double * c_mean = c_means[k];
        for(int d=0;d<D;d++)  c_mean[d] = c_tmpvec[d] / Nk;
      }
    }
  
  }
  

  bool MOG_diag_kmeans::dezombify_means() {

    static int counter = 0;

    bool zombie_mean = false;

    int k = 0;
    int max_count = count[k];
    int k_hog = k;
  
    for(int k=1;k<K;k++)  if( c_count[k] > max_count ) { max_count = c_count[k]; k_hog = k; }

    for(int k=0;k<K;k++) {
      if( c_count[k] == 0 ) {
      
        zombie_mean = true;
        if(verbose)  it_warning("MOG_diag_kmeans::dezombify_means(): detected zombie mean");
      
        if(k_hog == k) {
          it_warning("MOG_diag_kmeans::dezombify_means(): weirdness: k_hog == k");
          return(false);
        }
      
        if( counter >= c_count[k_hog] )  counter = 0;
      
        double * c_mean = c_means[k];
        double * c_x = c_X[ c_partitions[k_hog][counter] ];
      
        for(int d=0;d<D;d++)  c_mean[d] = 0.5*( c_means[k_hog][d] + c_x[d] );
        counter++;
      }

    }
    
    if(zombie_mean)  assign_to_means();

    return(true);
  }


  double MOG_diag_kmeans::measure_change() const {
  
    double tmp_dist = 0.0;
    for(int k=0;k<K;k++)   tmp_dist += dist( c_means[k], c_means_old[k] );
    return(tmp_dist);
  }


  void MOG_diag_kmeans::initial_means() {
  
    for(int d=0;d<D;d++) c_tmpvec[d] = 0.0;
  
    for(int n=0;n<N;n++) {
      double * c_x = c_X[n];
      for(int d=0;d<D;d++)  c_tmpvec[d] += c_x[d]; 
    }  
  
    for(int d=0;d<D;d++) c_tmpvec[d] /= N;
  
    int step = int(floor(double(N)/double(K)));
    for(int k=0;k<K;k++) {
      double * c_mean = c_means[k];
      double * c_x = c_X[k*step];
    
      for(int d=0;d<D;d++)  c_mean[d] = 0.5*(c_tmpvec[d] + c_x[d]);  
    }
  }


  void MOG_diag_kmeans::iterate() {
  
    for(int k=0;k<K;k++)  for(int d=0;d<D;d++)  c_means_old[k][d] = c_means[k][d];

    for(int i=0;i<max_iter;i++) {

      assign_to_means();
      if(!dezombify_means()) return;
      recalculate_means();

      double change = measure_change();

      if(verbose) std::cout << "MOG_diag_kmeans::iterate(): iteration = " << i << "  change = " << change << std::endl;
      if(change == 0) break;
  
      for(int k=0;k<K;k++)  for(int d=0;d<D;d++)   c_means_old[k][d] = c_means[k][d];
    }
  
  }
   
   
  void MOG_diag_kmeans::calc_means() {
    initial_means();
    iterate();
  }
   

  void MOG_diag_kmeans::calc_covs() {

    for(int k=0;k<K;k++) {
      int Nk = c_count[k];

      if(Nk >= 2) {
        double * c_mean = c_means[k];

        for(int d=0;d<D;d++) c_tmpvec[d] = 0.0;

        for(int n=0;n<Nk;n++) {
          double * c_x = c_X[ c_partitions[k][n] ];
          for(int d=0;d<D;d++) { double tmp = c_x[d] - c_mean[d];  c_tmpvec[d] += tmp*tmp; }
        }
      
        for(int d=0;d<D;d++)  c_diag_covs[k][d] = trust*(c_tmpvec[d] / (Nk-1.0) ) + (1.0-trust)*(1.0); 
      }
      else {
        for(int d=0;d<D;d++)  c_diag_covs[k][d] = 1.0;
      }
    }

  }


  void MOG_diag_kmeans::calc_weights() {
    for(int k=0;k<K;k++)  c_weights[k] = trust*(c_count[k] / double(N)) + (1.0-trust)*(1.0/K); 
    }


  void MOG_diag_kmeans::normalise_vectors() {

    for(int d=0;d<D;d++) {
      double acc = 0.0;  for(int n=0;n<N;n++)  acc += c_X[n][d];
      c_norm_mu[d] = acc / N;
    }

    for(int d=0;d<D;d++) {
      double acc = 0.0;  for(int n=0;n<N;n++) { double tmp = c_X[n][d] - c_norm_mu[d];  acc += tmp*tmp; } 
      c_norm_sd[d] = std::sqrt(acc/(N-1));
    }
  
    for(int n=0;n<N;n++) for(int d=0;d<D;d++) {
      c_X[n][d] -= c_norm_mu[d];
      if(c_norm_sd[d] > 0.0)  c_X[n][d] /= c_norm_sd[d];
    }
  }


  void MOG_diag_kmeans::unnormalise_vectors() {
  
    for(int n=0;n<N;n++)  for(int d=0;d<D;d++) {
      if(c_norm_sd[d] > 0.0)  c_X[n][d] *= c_norm_sd[d];
      c_X[n][d] += c_norm_mu[d];
    }
  }


  void MOG_diag_kmeans::unnormalise_means() {
  
    for(int k=0;k<K;k++) for(int d=0;d<D;d++) {
      if(norm_sd[d] > 0.0) c_means[k][d] *= c_norm_sd[d];
      c_means[k][d] += norm_mu[d];
    }
  }


  void MOG_diag_kmeans::run(MOG_diag &_model, Array<vec> &_X, int _max_iter=10, double _trust=0.5, bool _normalise=true) {
  
    it_assert( _model.is_valid(), "MOG_diag_kmeans::run(): given model is not valid" );
    it_assert( (_max_iter > 0), "MOG_diag_kmeans::run(): _max_iter needs to be greater than zero" );
    it_assert( ((_trust >= 0.0) && (_trust <= 1.0) ), "MOG_diag_kmeans::run(): _trust must be between 0 and 1 (inclusive)" );
  
    Array<vec> _means = _model.get_means(); Array<vec> _diag_covs = _model.get_diag_covs(); vec _weights = _model.get_weights();
    init( _means, _diag_covs, _weights );
    _means.set_size(0); _diag_covs.set_size(0); _weights.set_size(0);
  
    it_assert(check_size(_X), "MOG_diag_kmeans::run(): X is empty or contains vectors of wrong dimensionality" );
  
    N = _X.size();
  
    if(K > N)    it_warning("MOG_diag_kmeans::run(): K > N");
    else
    if(K > N/10) it_warning("MOG_diag_kmeans::run(): K > N/10");

    max_iter = _max_iter;
    trust = _trust;
  
    means_old.set_size(K);  for(int k=0;k<K;k++) means_old(k).set_size(D);
    partitions.set_size(K);  for(int k=0;k<K;k++) partitions(k).set_size(N);
    count.set_size(K);
    tmpvec.set_size(D);
    norm_mu.set_size(D);
    norm_sd.set_size(D);
  
    c_X = enable_c_access(_X);
    c_means_old = enable_c_access(means_old);
    c_partitions = enable_c_access(partitions);
    c_count = enable_c_access(count);
    c_tmpvec = enable_c_access(tmpvec);
    c_norm_mu = enable_c_access(norm_mu);
    c_norm_sd = enable_c_access(norm_sd);
      
    if(_normalise)  normalise_vectors();
  
    calc_means();  if(_normalise) { unnormalise_vectors(); unnormalise_means(); }  
    calc_covs();
    calc_weights();
  
    _model.init(means, diag_covs, weights);
  
    disable_c_access(c_X);
    disable_c_access(c_means_old);
    disable_c_access(c_partitions);
    disable_c_access(c_count);
    disable_c_access(c_tmpvec);
    disable_c_access(c_norm_mu);
    disable_c_access(c_norm_sd);
  
    means_old.set_size(0);
    partitions.set_size(0);
    count.set_size(0);
    tmpvec.set_size(0);
    norm_mu.set_size(0);
    norm_sd.set_size(0); 
  
    cleanup();
  
  }

}

