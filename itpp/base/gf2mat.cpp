/*!
 * \file
 * \brief Implementation of a class for algebra on GF(2) (binary) matrices
 * \author Erik G. Larsson
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

#include <itpp/base/gf2mat.h>
#include <itpp/base/specmat.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/converters.h>
#include <iostream>

namespace itpp {

  GF2mat::GF2mat(int i, int j) 
  {
    int k = j/(1<<lImax)+1;
    nrows=i;
    ncols=j;
    nwords=k;
    data.set_size(nrows,nwords);
    for (int i=0; i<nrows; i++) {
      for (int j=0; j<nwords; j++) {
	data(i,j) = 0;
      }
    }
  }

  GF2mat::GF2mat()   
  {
    nrows=1;
    ncols=1;
    nwords=1;
    data.set_size(1,1);
    data(0,0) = 0;
  }

  GF2mat::GF2mat(const bvec &x, bool row_or_column)
  {
    if (row_or_column==1) {     // create column vector
      nrows=length(x);
      ncols=1;
      nwords=1;
      data.set_size(nrows,1);
      for (int i=0; i<nrows; i++) {
	data(i,0) = 0;
	set(i,0,x(i));
      }
    } else {    // create row vector
      nrows=1;
      ncols=length(x);
      nwords=ncols/(1<<lImax)+1;
      data.set_size(1,nwords);
      for (int i=0; i<nwords; i++) {
	data(0,i) = 0;
      }
      for (int i=0; i<ncols; i++) {
	set(0,i,x(i));
      }
    }
  }
  

  GF2mat::GF2mat(const bmat &X)
  {
    it_error("not yet implemented");
  }


  GF2mat::GF2mat(const GF2mat_sparse &X)
  {
    nrows=X.rows();
    ncols=X.cols();
    //  GF2mat(m,n);
    int k = ncols/(1<<lImax)+1;
    nwords=k;
    data.set_size(nrows,nwords);
    for (int i=0; i<nrows; i++) {
      for (int j=0; j<nwords; j++) {
	data(i,j) = 0;
      }
    }

    for (int j=0; j<ncols; j++) {
      for (int i=0; i<X.get_col(j).nnz(); i++)   {
	bin b = X.get_col(j).get_nz_data(i);
	set(X.get_col(j).get_nz_index(i),j,b);
      }
    }
  }

  GF2mat::GF2mat(const GF2mat_sparse &X, int m1, int n1, int m2, int n2)
  {
    it_assert1(X.rows()>m2,"GF2mat()");
    it_assert1(X.cols()>n2,"GF2mat()");
    it_assert1(m1>=0 && n1>=0 && m2>=m1 && n2>=n1,"GF2mat::GF2mat()");
  
    nrows=m2-m1+1;
    ncols=n2-n1+1;
    nwords=ncols/(1<<lImax)+1;
    data.set_size(nrows,nwords);

    for (int i=0; i<nrows; i++) {
      for (int j=0; j<nwords; j++) {
	data(i,j) = 0;
      }
    }
  
    for (int i=0; i<nrows; i++) {
      for (int j=0; j<ncols; j++)   {
	bin b = X(i+m1,j+n1);
	set(i,j,b);
      }
    }
  }


  GF2mat::GF2mat(const GF2mat_sparse &X, ivec &columns)
  {
    it_assert1(X.cols()>max(columns),"GF2mat()");
    it_assert1(min(columns)>=0,"GF2mat()");
  
    nrows=X.rows();
    ncols=length(columns);
    nwords=ncols/(1<<lImax)+1;
    data.set_size(nrows,nwords);

    for (int i=0; i<nrows; i++) {
      for (int j=0; j<nwords; j++) {
	data(i,j) = 0;
      }
    }
  
    for (int j=0; j<ncols; j++) {
      for (int i=0; i<X.get_col(columns(j)).nnz(); i++)   {
	bin b = X.get_col(columns(j)).get_nz_data(i);
	set(X.get_col(columns(j)).get_nz_index(i),j,b);
      }
    }
  }

  GF2mat_sparse GF2mat::sparsify()
  {
    GF2mat_sparse Z(nrows,ncols);
    for (int i=0; i<nrows; i++) {
      for (int j=0; j<ncols; j++) {
	if (get(i,j)==1) {
	  Z.set(i,j,1);
	}
      }
    }
  
    return Z;
  }

  bvec GF2mat::bvecify()
  {
    it_assert1(nrows==1 || ncols==1,"GF2mat::bvecify()");
    int n=(nrows==1 ? ncols : nrows);
    bvec result(n);
    for (int i=0; i<n; i++) {
      result(i)=(nrows==1 ? get(0,i) : get(i,0));
    }
    return result;
  }


  void GF2mat::set_row(int i, bvec x)
  {
    it_assert1(length(x)==ncols,"GF2mat::set_row()");
    for (int j=0; j<ncols; j++) {
      set(i,j,x(j));
    }
  }
 
  void GF2mat::set_col(int j, bvec x)
  {
    it_assert1(length(x)==nrows,"GF2mat::set_col()");
    for (int i=0; i<nrows; i++) {
      set(i,j,x(i));
    }
  }


  GF2mat gf2dense_eye(int m)
  {
    GF2mat Z(m,m);
    for (int i=0; i<m; i++) {
      Z.set(i,i,1);
    }
    return Z;
  }

  GF2mat GF2mat::get_submatrix(int m1, int n1, int m2, int n2)
  {
    it_assert1(m1>=0 && n1>=0 && m2>=m1 && n2>=n1 && m2<nrows && n2<ncols,"GF2mat::submatrix");
    GF2mat result(m2-m1+1,n2-n1+1);

    for (int i=m1; i<=m2; i++) {
      for (int j=n1; j<=n2; j++) {
	result.set(i-m1,j-n1,get(i,j));
      }
    }

    return result;
  }


  GF2mat GF2mat::concatenate_vertical(const GF2mat &X)
  {
    it_assert1(X.ncols==ncols,"GF2mat::concatenate_vertical()");
    it_assert1(X.nwords==nwords,"GF2mat::concatenate_vertical()");
  
    GF2mat result(nrows+X.nrows,ncols);
    for (int i=0; i<nrows; i++) {
      for (int j=0; j<nwords; j++) {
	result.data(i,j) = data(i,j);
      }
    }

    for (int i=0; i<X.nrows; i++) {
      for (int j=0; j<nwords; j++) {
	result.data(i+nrows,j) = X.data(i,j);
      }
    }
  
    return result;
  }

  GF2mat GF2mat::concatenate_horizontal(const GF2mat &X)
  {
    //  std::cout << X.nrows << "     " << X.ncols << std::endl;
    //  std::cout << nrows << "     " << ncols << std::endl;
    it_assert1(X.nrows==nrows,"GF2mat::concatenate_horizontal()");

    GF2mat result(nrows,X.ncols+ncols);
    for (int i=0; i<nrows; i++) {
      for (int j=0; j<ncols; j++) {
	result.set(i,j,get(i,j));
      }
    }

    for (int i=0; i<nrows; i++) {
      for (int j=0; j<X.ncols; j++) {
	result.set(i,j+ncols,X.get(i,j));
      }
    }

    return result;
  }

  bvec GF2mat::get_row(int i)
  {
    bvec result = zeros_b(ncols);
    for (int j=0; j<ncols; j++) {
      result(j) = get(i,j);
    }

    return result;
  }

  bvec GF2mat::get_col(int j)
  {
    bvec result = zeros_b(nrows);
    for (int i=0; i<nrows; i++) {
      result(i) = get(i,j);
    }
  
    return result;
  }
  

  int GF2mat::T_fact(GF2mat &T, GF2mat &U, ivec &perm)
  {
    T = gf2dense_eye(nrows);
    U = *this;

    perm = zeros_i(ncols);
    for (int i=0; i<ncols; i++) {
      perm(i) = i;
    }
  
    if (nrows>250) {     // avoid cluttering output with this for small matrices...
      std::cout << "Performing T-factorization of GF(2) matrix...  rows: " << nrows << " cols: " 
	   << ncols << " .... " << std::endl;
      std::cout.flush();
    }
    int pdone=0;
    for (int j=0; j<nrows; j++) {
      // Now working on diagonal element j,j
      // First try find a row with a 1 in column i
      int i1,j1;
      for (i1=j; i1<nrows; i1++) {
	for (j1=j; j1<ncols; j1++) {
	  if (U.get(i1,j1)==1) { goto found; }
	}
      }

      return j;
    
    found: 
      U.swap_rows(i1,j);
      T.swap_rows(i1,j);
      U.swap_cols(j1,j);
    
      int temp = perm(j);
      perm(j) = perm(j1);
      perm(j1) = temp;
    
      // now subtract row i from remaining rows
      for (int i1=j+1; i1<nrows; i1++)  {
	if (U.get(i1,j)==1) { 
	  int ptemp = floor_i(100.0*(i1+j*nrows)/(nrows*nrows));
	  if (nrows>250 && ptemp>pdone+10) {     // this is only disturbing for small matrices
	    std::cout << ptemp << "% done." << std::endl;
	    std::cout.flush();
	    pdone=ptemp;
	  }
	  U.add_rows(i1,j);
	  T.add_rows(i1,j);
	}
      }
    }
    return nrows;
  }


  int GF2mat::T_fact_update_bitflip(GF2mat &T, GF2mat &U, ivec &perm, int rank, int r, int c)
  {
    // First, update U (before re-triangulization)
    int c0=c;
    for (c=0; c<ncols; c++) {
      if (perm(c)==c0) {
	goto foundidx;
      }
    }
    it_error("GF2mat::T_fact_update_bitflip() - internal error");  // should never reach this line

  foundidx:
    for (int i=0; i<nrows; i++) {
      if (T.get(i,r)==1) {
	U.addto_element(i,c,1);
      }
    }

    // first move column c to the end
    bvec lastcol = U.get_col(c);
    int temp_perm = perm(c);
    for (int j=c; j<ncols-1; j++) {
      U.set_col(j,U.get_col(j+1));
      perm(j) = perm(j+1);
    }
    U.set_col(ncols-1,lastcol);
    perm(ncols-1) = temp_perm;

    // then, if the matrix is tall, also move row c to the end 
    if (nrows>=ncols) {
      bvec lastrow_U = U.get_row(c);    
      bvec lastrow_T = T.get_row(c);
      for (int i=c; i<nrows-1; i++) {
	U.set_row(i,U.get_row(i+1));
	T.set_row(i,T.get_row(i+1));      
      }
      U.set_row(nrows-1,lastrow_U);
      T.set_row(nrows-1,lastrow_T);
    
      // Do Gaussian elimination on the last row 
      for (int j=c; j<ncols; j++)  {       
	if (U.get(nrows-1,j)==1) { 
	  U.add_rows(nrows-1,j);
	  T.add_rows(nrows-1,j);
	}
      }
    }
 
    // Now, continue T-factorization from the point (rank-1,rank-1)
    for (int j=rank-1; j<nrows; j++) { 
      int i1,j1;
      for (i1=j; i1<nrows; i1++) {
	for (j1=j; j1<ncols; j1++) {
	  if (U.get(i1,j1)==1) { 
	    goto found; 
	  }
	}
      }
    
      return j;
    
    found: 
      U.swap_rows(i1,j);
      T.swap_rows(i1,j);
      U.swap_cols(j1,j);
    
      int temp = perm(j);
      perm(j) = perm(j1);
      perm(j1) = temp;
    
      for (int i1=j+1; i1<nrows; i1++)  {
	if (U.get(i1,j)==1) { 
	  U.add_rows(i1,j);
	  T.add_rows(i1,j);
	}
      }
    }
  
    return nrows;
  }

  int GF2mat::T_fact_update_addcol(GF2mat &T, GF2mat &U, ivec &perm, bvec newcol)
  {
    int i0 = T.rows();
    int j0 = U.cols();
    it_assert1(j0>0,"GF2mat::T_fact_update_addcol()");
    it_assert1(i0==T.cols(),"GF2mat::T_fact_update_addcol()");
    it_assert1(i0==U.rows(),"GF2mat::T_fact_update_addcol()");
    it_assert1(length(perm)==j0,"GF2mat::T_fact_update_addcol()");
    it_assert1(U.get(j0-1,j0-1)==1,"GF2mat::T_fact_update_addcol()");
    it_assert0(U.row_rank()==j0,"GF2mat::T_fact_update_addcol()");  // VERY TIME CONSUMING TEST

    bvec z = T*newcol;
    GF2mat Utemp = U.concatenate_horizontal(GF2mat(z,1));

    // start working on position (j0,j0)
    int i;
    for (i=j0; i<i0; i++) {
      if (Utemp.get(i,j0)==1) {
	goto found;
      }
    }
    return 0; // adding the new column would not improve the rank

  found:
    perm.set_length(j0+1,true);
    perm(j0) = j0;

    Utemp.swap_rows(i,j0);
    T.swap_rows(i,j0);
  
    for (int i1=j0+1; i1<i0; i1++)  {
      if (Utemp.get(i1,j0)==1) { 
	Utemp.add_rows(i1,j0);
	T.add_rows(i1,j0);
      }
    }

    U = Utemp;
    return 1; // the new column was successfully added
  }




  GF2mat GF2mat::inverse()
  {
    it_assert1(nrows==ncols,"GF2mat::inverse(): Matrix must be square");

    // first compute the T-factorization
    GF2mat T,U;
    ivec perm;
    int rank = T_fact(T,U,perm);
    it_assert1(rank==ncols,"GF2mat::inverse(): Matrix is not full rank"); 
 
    // backward substitution
    for (int i=ncols-2; i>=0; i--) {
      for (int j=ncols-1; j>i; j--) {
	if (U.get(i,j)==1) {
	  U.add_rows(i,j);
	  T.add_rows(i,j);
	}
      }
    }
    T.permute_rows(perm,1);
    return T;
  }

  int GF2mat::row_rank()
  {
    GF2mat T,U;
    ivec perm;
    int rank = T_fact(T,U,perm);
    return rank;  
  }

  bool GF2mat::is_zero()
  {
    for (int i=0; i<nrows; i++) {
      for (int j=0; j<nwords; j++) {
	if (data(i,j)!=0) { 
	  return false; 
	}
      }
    }
    return true;
  }

  bool GF2mat::operator==(const GF2mat &X) const
  {
    if (X.nrows!=nrows) { return false; }
    if (X.ncols!=ncols) { return false; }
    it_assert1(X.nwords==nwords,"GF2mat::operator==()");
 
    for (int i=0; i<nrows; i++) {
      for (int j=0; j<nwords; j++) {
	//      if (X.get(i,j)!=get(i,j)) { 
	if (X.data(i,j)!=data(i,j)) {
	  return false; 
	}
      }
    }
    return true;  
  }


  void GF2mat::add_rows(int i, int j)
  {
    it_assert0(i>=0 && i<nrows,"GF2mat::add_rows() - out of range");
    it_assert0(j>=0 && j<nrows,"GF2mat::add_rows() - out of range");
    for (int k=0; k<nwords; k++)   {
      data(i,k) ^= data(j,k);    
    }  
  }

  void GF2mat::swap_rows(int i, int j)
  {
    it_assert0(i>=0 && i<nrows,"GF2mat::swap_rows()");
    it_assert0(j>=0 && j<nrows,"GF2mat::swap_rows()");
    for (int k=0; k<nwords; k++) {
      int temp = data(i,k);
      data(i,k) = data(j,k);
      data(j,k) = temp;
    }
  }

  void GF2mat::swap_cols(int i, int j)
  {
    it_assert0(i>=0 && i<ncols,"GF2mat::swap_cols()");
    it_assert0(j>=0 && j<ncols,"GF2mat::swap_cols()");
    bvec temp = get_col(i);
    set_col(i,get_col(j));
    set_col(j,temp);
  }


  void GF2mat::operator=(const GF2mat &X)
  {
    nrows=X.nrows;
    ncols=X.ncols;
    nwords=X.nwords;
    data = X.data;
  }

  GF2mat operator*(const GF2mat &X, const GF2mat &Y)
  {
    it_assert1(X.ncols==Y.nrows,"GF2mat::operator*");
    it_assert0(X.nwords>0,"Gfmat::operator*");
    it_assert0(Y.nwords>0,"Gfmat::operator*");

    /*  
    // this can be done more efficiently?
    GF2mat result(X.nrows,Y.ncols);
    for (int i=0; i<X.nrows; i++) {
    for (int j=0; j<Y.ncols; j++) {
    bin b=0;
    for (int k=0; k<X.ncols; k++) {
    bin x = X.get(i,k);
    bin y = Y.get(k,j);
    b ^= (x&y);
    }
    result.set(i,j,b);
    }
    }
    return result; */

    // is this better?
    return mult_trans(X,Y.transpose());
  }

  bvec operator*(const GF2mat &X, const bvec &y) 
  {
    it_assert1(length(y)==X.ncols,"GF2mat::operator*");
    it_assert0(X.nwords>0,"Gfmat::operator*");

    /*  
    // this can be done more efficiently?
    bvec result = zeros_b(X.nrows);
    for (int i=0; i<X.nrows; i++) {
    // do the nwords-1 data columns first
    for (int j=0; j<X.nwords-1; j++) {
    int ind = j<<lImax;
    unsigned short r=X.data(i,j); // NB. THIS MUST BE UNSIGNED FOR >> TO BE WELL DEFINED
    while (r) {
    result(i) ^= (r & y(ind));
    r >>= 1;    
    ind++;
    }
    }
    // do the last column separately
    for (int j=(X.nwords-1)<<lImax; j<X.ncols; j++) {
    result(i) ^= (X.get(i,j) & y(j));
    }
    }      
    return result; */

    // is this better?
    return (mult_trans(X,GF2mat(y,0))).bvecify();
  }

  GF2mat mult_trans(const GF2mat &X, const GF2mat &Y)
  {
    it_assert1(X.ncols==Y.ncols,"GF2mat mult_trans");
    it_assert0(X.nwords>0,"GF2mat mult_trans");
    it_assert0(Y.nwords>0,"GF2mat mult_trans");
    it_assert0(X.nwords==Y.nwords,"GF2mat mult_trans");
  
    GF2mat result(X.nrows,Y.nrows);

    for (int i=0; i<X.nrows; i++) {
      for (int j=0; j<Y.nrows; j++) {
	bin b=0;
	int kindx =i; 
	int kindy =j;
	for (int k=0; k<X.nwords; k++) {
	  //	unsigned short int r=(X.data(i,k) & Y.data(j,k));
	  unsigned short int r=X.data(kindx) & Y.data(kindy);
	  /* The following can be speeded up by using a small lookup
	     table for the number of ones and shift only a few times (or
	     not at all if table is large) */
	  while (r) {
	    b ^= r&1;
	    r>>=1;
	  };
	  kindx += X.nrows;
	  kindy += Y.nrows;
	}
	result.set(i,j,b);
      }
    }
    return result;
  }

  GF2mat GF2mat::transpose() const
  {
    // CAN BE SPEEDED UP
    GF2mat result(ncols,nrows);
  
    for (int i=0; i<nrows; i++) {
      for (int j=0; j<ncols; j++) {
	result.set(j,i,get(i,j));
      }
    }
    return result;
  }

  GF2mat operator+(const GF2mat &X, const GF2mat &Y)
  {
    it_assert1(X.nrows==Y.nrows,"GF2mat operator+");
    it_assert1(X.ncols==Y.ncols,"GF2mat operator+");
    it_assert1(X.nwords==Y.nwords,"GF2mat operator+");
    GF2mat result(X.nrows,X.ncols);
  
    for (int i=0; i<X.nrows; i++) {
      for (int j=0; j<X.nwords; j++) {
	result.data(i,j) = X.data(i,j) ^ Y.data(i,j);
      }
    }

    return result;
  }

  void GF2mat::permute_cols(ivec &perm, bool I)
  {
    it_assert1(length(perm)==ncols,"GF2mat::permute_cols()");

    GF2mat temp = (*this);
    for (int j=0; j<ncols; j++) {
      if (I==0) {
	set_col(j,temp.get_col(perm(j)));
      }  else { 
	set_col(perm(j),temp.get_col(j));
      }
    }
  }

  void GF2mat::permute_rows(ivec &perm, bool I)
  {
    it_assert1(length(perm)==nrows,"GF2mat::permute_rows()");

    GF2mat temp = (*this);
    for (int i=0; i<nrows; i++) {
      if (I==0) {
	for (int j=0; j<nwords; j++) {
	  data(i,j) = temp.data(perm(i),j);
	}
	//      set_row(i,temp.get_col(perm(i)));
      }  else { 
	for (int j=0; j<nwords; j++) {
	  data(perm(i),j) = temp.data(i,j);
	}
	//      set_row(perm(i),temp.get_row(i));
      }
    }
  }


  std::ostream &operator<<(std::ostream &os, const GF2mat &X)
  {
    int i,j;
    os << "---- GF(2) matrix of dimension " << X.nrows << "*" << X.ncols
       << " -- Density: " << X.density() << " ----" << std::endl;

    for (i=0; i<X.nrows; i++) {
      std::cout << "      ";
      for (j=0; j<X.ncols; j++) {
	os << X.get(i,j) << " ";
      }
      os << std::endl;
    }

    return os;
  }

  double GF2mat::density() const
  {
    int no_of_ones=0;

    for (int i=0; i<nrows; i++) {
      for (int j=0; j<ncols; j++) {
	no_of_ones += (get(i,j)==1 ? 1 : 0);
      }
    }

    return ((double) no_of_ones)/(nrows*ncols);
  }


  it_file &operator<<(it_file &f, const GF2mat &X)
  {
    int bytecount = 3*sizeof(int)+X.nrows*X.nwords*sizeof(int);
    f.write_data_header("GF2mat", bytecount);

    f.low_level_write(X.nrows);
    f.low_level_write(X.ncols);
    f.low_level_write(X.nwords);
    for (int i=0; i<X.nrows; i++) {
      for (int j=0; j<X.nwords; j++) {
	short r=X.data(i,j);
	f.low_level_write(r);
      }
    }
    return f; 
  }

  it_ifile &operator>>(it_ifile &f, GF2mat &X)
  {
    it_file::data_header h;
  
    f.read_data_header(h);
    if (h.type == "GF2mat") {
      f.low_level_read(X.nrows);
      f.low_level_read(X.ncols);
      f.low_level_read(X.nwords);
      X.data.set_size(X.nrows,X.nwords);
      for (int i=0; i<X.nrows; i++) {
	for (int j=0; j<X.nwords; j++) {
	  short r;
	  f.low_level_read(r);
	  X.data(i,j) = ((unsigned short int) r);
	}
      }
    
    }  else {
      //    throw it_file_base::Error("Wrong type");
      it_error("it_ifile &operator>>() - internal error");
    }

    return f;
  }

} // namespace itpp

