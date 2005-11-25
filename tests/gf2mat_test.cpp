/*!
 * \file
 * \brief Test program for a class for algebra on GF(2) (binary) matrices
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
 * Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
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

/*
  NOTE:

  - With the debugging options enables at compilation, this program
    takes a while to run.
*/

#include <itpp/itbase.h>

using std::cout;
using std::endl;
using namespace itpp;

GF2mat random_matrix(int m, int n)
{
  GF2mat Z(m,n);
  for (int j=0; j<n; j++) {
    for (int i=0; i<m; i++) {
      Z.set(i,j,(rand()%2==0 ? 1 : 0));
    } 
  }
  return Z;
};

int main()
{

  // Test of file I/O
  GF2mat Z = random_matrix(301,179);

  it_file f1("gf2mat_test.it");
  f1 << Name("Z") << Z;
  f1.close();
  cout << "Saved successfully." << endl;
  
  it_ifile f2("gf2mat_test.it");
  GF2mat Z_temp;
  f2 >> Name("Z") >> Z_temp;
  f2.close();
  it_assert(Z==Z_temp,"error");
  cout << "Read successfully." << endl;

  // Binary vector
  bvec b=randb(Z.cols());
  it_assert(GF2mat(Z*b,1)==Z*GF2mat(b,1),"error");
  cout << "Passed matrix-vector multiplication test" << endl;
  
  // Multiplication test
  GF2mat W = random_matrix(139,Z.rows());
  GF2mat temp1=W*Z;
  cout << "computed product." << endl;
  GF2mat temp2=GF2mat(W.sparsify()*Z.sparsify());
  it_assert(temp1==temp2,"error");
  Z=Z.transpose();
  it_assert(W*Z.transpose()==mult_trans(W,Z),"error");
  cout << "Passed matrix-matrix multiplication test" << endl;
  
  // Transpose
  it_assert(GF2mat(b,0)==GF2mat(b,1).transpose(),"error");
  it_assert(GF2mat(b,1)==GF2mat(b,0).transpose(),"error");
  GF2mat Y = random_matrix(Z.cols(),73);
  it_assert((Z*Y).transpose()==Y.transpose()*Z.transpose(),"error");
  cout << "Passed transpose test." << endl;

  // Concatenation
  int m=Z.rows();
  int n=Z.cols();
  it_assert(Z == Z.get_submatrix(0,0,m-1,27).concatenate_horizontal(Z.get_submatrix(0,28,m-1,n-1)),"error");
  it_assert(Z == Z.get_submatrix(0,0,13,n-1).concatenate_vertical(Z.get_submatrix(14,0,m-1,n-1)),"error");
  cout << "Passed concatenation test." << endl;
  
  // Assignment operator
  GF2mat P=Z;
  it_assert(P==Z,"error");
  it_assert((P+Z).is_zero(),"error");
  cout << "Passed assignment operator test." << endl;

  // Sparse-dense conversions
  GF2mat_sparse A(Z.rows(),Z.cols());
  for (int i=0; i<Z.rows(); i++) {
    for (int j=0; j<Z.cols(); j++) {
      if (Z.get(i,j)==1) { 
	A.set(i,j,1);
      }
    }
  }
  it_assert(GF2mat(A)==Z,"error");
  GF2mat_sparse C = Z.sparsify();
  it_assert(C.full()==A.full(),"error");
  cout << "Passed sparse test." << endl;

  Z = random_matrix(100,75);

  // Get rows and columns
  cout << "Z.get_row(1)=" << Z.get_row(1) << endl;
  cout << "Z.get_col(2)=" << Z.get_col(2) << endl;

  // Print a submatrix on the screen
  cout << "Z.get_submatrix(1,1,6,4): " << Z.get_submatrix(1,1,6,4) << endl;


  // Test of T-factorization
  int dim=250;
  
  for (int trial=0; trial<100; trial++) {
    cout << "Testing T-factorization, realization: " << trial << endl;
    GF2mat X=random_matrix(rand()%dim+1,rand()%dim+1);
    GF2mat T,U;
    ivec perm;
    X.T_fact(T,U,perm);
    GF2mat W = T*X;
    W.permute_cols(perm,0);
    it_assert(U==W,"error");
  }

  // Test of inversion
  for (int trial=0; trial<100; trial++) {
    cout << "Testing inversion, realization: " << trial << endl;
    GF2mat X=random_matrix(dim,dim);
    while (X.row_rank()!=dim) {
      X.set(rand()%dim,rand()%dim,rand()%2);
    }    
    it_assert(X*X.inverse() == gf2dense_eye(dim),"error");
    it_assert(X.inverse()*X == gf2dense_eye(dim),"error");
  }
  
  // Test of the T-factorization bitflip update
  for (int trial=0; trial<100; trial++) {
    cout << "Testing the T-factorization bitflip update, realization: " << trial;
    GF2mat X=random_matrix(rand()%dim+1,rand()%dim+1);
    cout << " dimension: " << X.rows() << "*" << X.cols();
    GF2mat T,U;
    ivec perm;
    int rank = X.T_fact(T,U,perm);
    cout << " rank:" << rank << endl;

    GF2mat Tnew=T;
    GF2mat Unew=U;
    ivec permnew = perm;
    for (int trial2=0; trial2<10; trial2++) {
      int i = rand()%X.rows();
      int j = rand()%X.cols();
      X.addto_element(i,j,1);
      X.T_fact_update_bitflip(Tnew,Unew,permnew,rank,i,j);
 
      GF2mat W = Tnew*X;
      W.permute_cols(permnew,0);
      it_assert(Unew==W,"error");
    }
  }


  // Test of the T-factorization add-column update
  for (int trial=0; trial<100; trial++) {
    cout << "Testing the T-factorization add-column update, realization number: " << trial << endl;
    bvec c = randb(dim);
    GF2mat X(c,1);

    GF2mat T,U;
    ivec perm;
    X.T_fact(T,U,perm);

    for (int trial2=0; trial2<100; trial2++) {
      bvec c = randb(dim);
      //      cout << X << endl;
      //      cout << GF2mat(c,1) << endl;
      GF2mat Xtemp = X.concatenate_horizontal(GF2mat(c,1));
      int success = Xtemp.T_fact_update_addcol(T,U,perm,c);
      if (success==1)	{ 
	X=Xtemp; 
      }
      //      cout << "rank was: " << X.row_rank() << endl;

      GF2mat W = T*X;
      W.permute_cols(perm,0);
      it_assert(U==W,"error");
    }
  }


  cout << "All tests successfully passed." << endl;
}
