/*!
 * \file 
 * \brief diagonal Mixture of Gaussians class - source file
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

#include <itpp/stat/mog_diag.h>


namespace itpp {

  void MOG_diag::zero_all_ptrs() {
    c_means             = 0;
    c_diag_covs         = 0;
    c_diag_covs_inv_etc = 0;
    c_weights           = 0;
    c_log_weights       = 0;
    c_log_det_etc       = 0;
    c_tmpvecK           = 0;
  }


  void MOG_diag::free_all_ptrs() {
    c_means             = disable_c_access(c_means);
    c_diag_covs         = disable_c_access(c_diag_covs);
    c_diag_covs_inv_etc = disable_c_access(c_diag_covs_inv_etc);
    c_weights           = disable_c_access(c_weights);
    c_log_weights       = disable_c_access(c_log_weights);
    c_log_det_etc       = disable_c_access(c_log_det_etc);
    c_tmpvecK           = disable_c_access(c_tmpvecK);
  }


  void MOG_diag::setup_means() {
    MOG_generic::setup_means();
    disable_c_access(c_means);
    c_means = enable_c_access(means);
  }


  void MOG_diag::setup_covs() {
    MOG_generic::setup_covs();
    if(full) return;
  
    disable_c_access(c_diag_covs);
    disable_c_access(c_diag_covs_inv_etc);
    disable_c_access(c_log_det_etc);
  
    c_diag_covs         = enable_c_access(diag_covs);
    c_diag_covs_inv_etc = enable_c_access(diag_covs_inv_etc);
    c_log_det_etc       = enable_c_access(log_det_etc);
  }


  void MOG_diag::setup_weights() {
    MOG_generic::setup_weights();
  
    disable_c_access(c_weights);
    disable_c_access(c_log_weights);
  
    c_weights = enable_c_access(weights);
    c_log_weights = enable_c_access(log_weights);
  }


  void MOG_diag::setup_misc() {
    disable_c_access(c_tmpvecK);
    tmpvecK.set_size(K);
    c_tmpvecK = enable_c_access(tmpvecK);
    
    MOG_generic::setup_misc();
    if(full) convert_to_diag_internal();
  }


  void MOG_diag::load(const std::string &_name) {
    MOG_generic::load(_name);
    if(full) convert_to_diag();
  }


  double ** MOG_diag::enable_c_access(Array<vec> & _A) {
    int rows = _A.size();
    double ** A = (double **)std::malloc(rows*sizeof(double *));
    if(A)  for(int row=0;row<rows;row++)  A[row] = _A(row)._data();
    return(A);
  }

  int ** MOG_diag::enable_c_access(Array<ivec> & _A) {
    int rows = _A.size();
    int ** A = (int **)std::malloc(rows*sizeof(int *));
    if(A)  for(int row=0;row<rows;row++)  A[row] = _A(row)._data();
    return(A);
  }

  double ** MOG_diag::disable_c_access(double ** A) { if(A) std::free(A); return(0); }
  int ** MOG_diag::disable_c_access(int ** A) { if(A) std::free(A); return(0); }

  double * MOG_diag::enable_c_access(vec & _v) { return _v._data(); }
  int * MOG_diag::enable_c_access(ivec & _v) { return _v._data(); }

  double * MOG_diag::disable_c_access(double * v) { return(0); }
  int * MOG_diag::disable_c_access(int * v) { return(0); }

}
