/*!
 * \file 
 * \brief Implementation of a Low-Density Parity Check (LDPC) codec.
 * \author Erik G. Larsson and Mattias Andersson
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

#include <itpp/comm/ldpc.h>

using std::cout;
using std::endl;

namespace itpp {

  // ---------------------------------------------------------------------------
  // LDPC_Parity_Matrix
  // ---------------------------------------------------------------------------
  
  template <class T>
  ivec LDPC_Parity_Matrix::index_nonzeros(Sparse_Vec<T> x) const
  {
    ivec r;
    for (int i=0; i<x.nnz(); i++)     {
      if (x.get_nz_data(i)!=0) 	{
	r = concat(r,x.get_nz_index(i));
      }
    }
    return r;
  }

  void LDPC_Parity_Matrix::initialize(int nc, int nv)
  {
    ncheck=nc;
    nvar=nv;
    H=GF2mat_sparse(ncheck,nvar);
    Ht=GF2mat_sparse(nvar,ncheck);
    sumX1 = zeros_i(nvar);
    sumX2 = zeros_i(ncheck);
  }

  LDPC_Parity_Matrix::LDPC_Parity_Matrix()
  {
    initialize(1,1);
  }

  LDPC_Parity_Matrix::LDPC_Parity_Matrix(int nc, int nv) 
  {
    initialize(nc,nv);
  }

  LDPC_Parity_Matrix::LDPC_Parity_Matrix(std::string filename, std::string format)
  {
    if (format=="alist") {
      *this = LDPC_Parity_Matrix(GF2mat_sparse_alist(filename));
    } else { 
      it_error("not implemented");
    };
  }

  LDPC_Parity_Matrix::LDPC_Parity_Matrix(const GF2mat_sparse_alist &alist) 
  {
    GF2mat_sparse X=alist.to_sparse();

    initialize(X.rows(),X.cols());
    // brute force copy from X to this
    for (int i=0; i<ncheck; i++) {
      for (int j=0; j<nvar; j++) {
 	if (X(i,j)) {
 	  set(i,j,1);
 	}
      }
    }    
  }

  void LDPC_Parity_Matrix::set(int i, int j, bin x) 
  {
    it_assert0((x==0) | (x==1),"LDPC_Parity_Matrix::set");
    it_assert0((i>=0) & (i<ncheck),"LDPC_Parity_Matrix::set");
    it_assert0((j>=0) & (j<nvar),"LDPC_Parity_Matrix::set");
    it_assert0(H(i,j)==Ht(j,i),"LDPC_Parity_Matrix internal error");

    int old_element = H(i,j);
    if (x==1)   {
      H.set(i,j,1);
      Ht.set(j,i,1);
    }
    else    {
      H.clear_elem(i,j);
      Ht.clear_elem(j,i);
    }
    int diff = ((int) x)-((int) old_element);
    sumX1(j) += diff;
    sumX2(i) += diff;
  
    it_assert0(H(i,j)==x,"LDPC_Parity_Matrix::set");
    it_assert0(Ht(j,i)==x,"LDPC_Parity_Matrix::set");
  }

  //   ivec LDPC_Parity_Matrix::get_rowdegree()  const
  //   {
  //     ivec rdeg = zeros_i(Nmax);
  //     for (int i=0; i<ncheck; i++)     {
  //       rdeg(sumX2(i))++;
  //     }  
  //     return rdeg;
  //   };
  
  //   ivec LDPC_Parity_Matrix::get_coldegree()  const
  //   {
  //     ivec cdeg = zeros_i(Nmax);
  //     for (int j=0; j<nvar; j++)     {
  //       cdeg(sumX1(j))++;
  //     }
  //     return cdeg;
  //   };
  
  void LDPC_Parity_Matrix::display_stats() const
  {
    int cmax=max(sumX1);
    int vmax=max(sumX2);
    vec vdeg = zeros(cmax+1); //no of var nodes with n neighbours
    vec cdeg = zeros(vmax+1); //no of check nodes with n neighb
    for (int col=0; col<nvar; col++){
      vdeg(length(index_nonzeros(get_col(col))))++;  
      it_assert(sumX1(col)==length(index_nonzeros(get_col(col))),"internal error");
    }
    for (int row=0; row<ncheck; row++){
      cdeg(length(index_nonzeros(get_row(row))))++;
      it_assert(sumX2(row)==length(index_nonzeros(get_row(row))),"internal error");
    }
    
    //from edge perspective
    //number of edges connected to vnodes of degree n
    vec vdegedge = elem_mult(vdeg,linspace(0,vdeg.length()-1,
					   vdeg.length()));
    //-"- but cnodes 
    vec cdegedge = elem_mult(cdeg,linspace(0,cdeg.length()-1,
					   cdeg.length()));

    int edges = sum(elem_mult(to_ivec(linspace(0,vdeg.length()-1,
					       vdeg.length())),
			      to_ivec(vdeg)));
  
    cout << "--- LDPC parity check matrix ---" << endl;
    cout << "Dimension: ncheck*nvar = " << ncheck << " * " << nvar << endl;
    cout << "Variable node degree distribution from node perspective:" 
	 << endl << vdeg/nvar << endl;
    cout << "Check node degree distribution from node perspective:" 
	 << endl << cdeg/ncheck << endl;
    cout << "Variable node degree distribution from edge perspective:" 
	 << endl << vdegedge/edges << endl;
    cout << "Check node degree distribution from edge perspective:" 
	 << endl << cdegedge/edges << endl;  
    cout << "--------------------------------" << endl;
  }


  GF2mat_sparse_alist LDPC_Parity_Matrix::export_alist() const 
  {
    GF2mat_sparse_alist alist;
    alist.from_sparse(H);
    return alist;
  }

  int LDPC_Parity_Matrix::cycle_removal_MGW(int Maxcyc)
  {
    typedef Sparse_Mat<short> Ssmat;
    typedef Sparse_Vec<short> Ssvec;
    
    Maxcyc -= 2; 

    // Construct the adjacency matrix of the graph
    Ssmat G(ncheck+nvar,ncheck+nvar,5);
    for (int j=0; j<nvar; j++) {
      GF2vec_sparse col = get_col(j);
      for (int i=0; i<col.nnz(); i++) {
	if (get(col.get_nz_index(i),j)==1) {
	  G.set(col.get_nz_index(i),j+ncheck,1);
	  G.set(ncheck+j,col.get_nz_index(i),1);
	}
      }
    }

    Array<Ssmat> Gpow(Maxcyc);  
    Gpow(0).set_size(ncheck+nvar,ncheck+nvar,1);
    Gpow(0).clear();
    for (int i=0; i<ncheck+nvar; i++) {
      Gpow(0).set(i,i,1); 
    }
    Gpow(1) = G;

    /* Main cycle elimination loop starts here. Note that G and all
       powers of G are symmetric matrices. This fact is exploited in
       the code.
    */
    int r;
    int cycles_found =0;
    int scl=Maxcyc;
    for (r=4; r<=Maxcyc; r+=2) {
      // compute the next power of the adjacency matrix
      Gpow(r/2) = Gpow(r/2-1)*G;            
      bool traverse_again;
      do {
	traverse_again=false;
	cycles_found=0;
	cout << "Starting new pass of cycle elimination, target girth "
	     << (r+2) << "... " << endl;
	int pdone=0;
	for (int j=0; j<ncheck+nvar; j++) { // loop over elements of G
	  for (int i=0; i<ncheck+nvar; i++ ) {
	    int ptemp = floor_i(100.0*(i+j*(ncheck+nvar))/
				((nvar+ncheck)*(nvar+ncheck)));
	    if (ptemp>pdone+10) {
	      cout << ptemp << "% done." << endl;
	      cout.flush();
	      pdone=ptemp;
	    }
	  
	    if (((Gpow(r/2))(i,j) >= 2)  && ((Gpow(r/2-2))(i,j)==0)) {
	      // Found a cycle.
	      cycles_found++;

	      // choose k
	      ivec tmpi = index_nonzeros(elem_mult(Gpow(r/2-1).get_col(i),
						   G.get_col(j)));
// 	      int k = tmpi(rand()%length(tmpi));    
	      int k = tmpi(randi(0,length(tmpi)-1));    
	      it_assert0(G(j,k)==1 && G(k,j)==1,
			 "LDPC_Parity_Matrix::remove_cycles_MGW() internal error");

	      // determine candidate edges for an edge swap
	      Ssvec rowjk = Gpow(r/2)*(Gpow(r/2-1).get_col(j)
				       +Gpow(r/2-1).get_col(k));
	      int p,l;
	      ivec Ce_ind = sort_index(randu(nvar+ncheck)); // random order

	      for (int s=0; s<nvar+ncheck; s++)  {
		l = Ce_ind(s);
		if (rowjk(l)!=0) { continue; }
		ivec colcandi = index_nonzeros(G.get_col(l));
		if (length(colcandi)>0)  {
		  // select a node p which is a member of Ce
		  for (int u=0; u<length(colcandi); u++) {
		    p = colcandi(u);
		    if (p!=l) {
		      if (rowjk(p)==0) {
			goto found_candidate_vector;
		      }
		    }
		  }
		}
	      }
	      continue; // go to the next entry (i,j)
	      
	    found_candidate_vector:
	      // swap edges

	      if (p>=ncheck) { int z=l; l=p; p=z; } 
	      if (j>=ncheck) { int z=k; k=j; j=z; } 

	      // Swap endpoints of edges (p,l) and (j,k)
	      // cout << "(" << j << "," << k << ")<->(" 
	      // << p << "," << l << ") " ;
	      // cout << ".";
	      // cout.flush();

	      // Update the matrix
	      it_assert0((get(j,k-ncheck)==1) && (get(p,l-ncheck)==1),
			 "LDPC_Parity_Matrix::remove_cycles_MGW() internal error");
	      set(j,k-ncheck,0);
	      set(p,l-ncheck,0);
	      it_assert0((get(j,l-ncheck)==0) && (get(p,k-ncheck)==0),
			 "LDPC_Parity_Matrix::remove_cycles_MGW() internal error");
	      set(j,l-ncheck,1);
	      set(p,k-ncheck,1);
	      
	      // cout << "Number of cycles remaining: " 
	      // << check_for_cycles(r) << endl;

	      // Update adjacency matrix
	      it_assert0(G(p,l)==1 && G(l,p)==1 && G(j,k)==1 
			 && G(k,j)==1,"G");
	      it_assert0(G(j,l)==0 && G(l,j)==0 && G(p,k)==0 
			 && G(k,p)==0,"G");
      	      
	      // Delta is the update term to G
	      Ssmat Delta(ncheck+nvar,ncheck+nvar,2);    
	      Delta.set(j,k,-1);    Delta.set(k,j,-1);
	      Delta.set(p,l,-1);    Delta.set(l,p,-1);     
	      Delta.set(j,l,1);	    Delta.set(l,j,1);
	      Delta.set(p,k,1);	    Delta.set(k,p,1); 

	      // update G and its powers
	      G = G+Delta;
	      it_assert0(G(p,l)==0 && G(l,p)==0 && G(j,k)==0 
			 && G(k,j)==0,"G");
	      it_assert0(G(j,l)==1 && G(l,j)==1 && G(p,k)==1 
			 && G(k,p)==1,"G");

	      Gpow(1)=G;	
	      Gpow(2)=G*G;
	      for (int z=3; z<=r/2; z++) {
		Gpow(z) = Gpow(z-1)*G;
		//	      cout << "power: " << z << "    density: " 
		// << Gpow(z).density() << endl;
	      }

	      traverse_again=true;
	    } // if G()...
	  } // loop over i
	}  // loop over j
	if ((!traverse_again) && (cycles_found>0)) { // no point continue
	  scl=r-2;
	  goto finished;
	}
      }  while (cycles_found!=0);
      scl=r;  // there were no cycles of length r; move on to next r
      cout << "Achieved girth " << (scl+2) 
	   << ". Proceeding to next level." << endl;
    } // loop over r

  finished:
    int girth=scl+2;  // scl=length of smallest cycle
    cout << "Cycle removal (MGW algoritm) finished. Graph girth: " 
	 << girth << ". Cycles remaining on next girth level: " 
	 << cycles_found << endl;
 
    return girth;
  } 

  void LDPC_Parity_Matrix::generate_regular_ldpc(int Nvar, int k, int l, 
						 std::string method,
						 ivec options)
  {
    int Ncheck_actual = static_cast<int>( round( static_cast<double>(Nvar) *
						 static_cast<double>(k) /
						 static_cast<double>(l) ) );
    // C, R: Target number of columns/rows with certain number of ones
    ivec C = zeros_i(Nmax);
    ivec R = zeros_i(Nmax);
    C(k) = Nvar;
    R(l) = Ncheck_actual;

    // ---------------

    if (method=="rand") {
      generate_random_H(C,R,options);
    } else { 
      it_error("not implemented"); 
    };
  }

  void LDPC_Parity_Matrix::generate_irregular_ldpc(int Nvar, vec var_deg, 
						   vec chk_deg, std::string method,
						   ivec options)
  {
    // compute the degree distributions from a node perspective
    vec Vi = linspace(1,length(var_deg),length(var_deg)); 
    vec Ci = linspace(1,length(chk_deg),length(chk_deg)); 
    // Compute number of cols with n 1's
    // C, R: Target number of columns/rows with certain number of ones
    ivec C = to_ivec(round(Nvar*elem_div(var_deg,Vi)
			   /sum(elem_div(var_deg,Vi))));   
    C = concat(0,C);
    int edges = sum(elem_mult(to_ivec(linspace(0,C.length()-1,
					       C.length())),C));
    ivec R = to_ivec(round(edges*elem_div(chk_deg,Ci))); 
    R = concat(0,R);
    vec Ri = linspace(0,length(R)-1,length(R));
    vec Coli = linspace(0,length(C)-1,length(C));

    // trim to deal with inconsistencies due to rounding errors
    if (sum(C)!=Nvar) {
      ivec ind = find(C==max(C));
      C(ind(0)) = C(ind(0)) - (sum(C)-Nvar);
    }

    //the number of edges calculated from R must match the number of
    //edges calculated from C
    while (sum(elem_mult(to_vec(R),Ri)) != 
	   sum(elem_mult(to_vec(C),Coli))) {
      //we're only changing R, this is probably(?) better for irac codes
      if (sum(elem_mult(to_vec(R),Ri)) > sum(elem_mult(to_vec(C),Coli))) {
	//remove an edge from R
	ivec ind = find(R == max(R));
	int old = R(ind(0));
	R.set(ind(0),old-1);
	old = R(ind(0)-1);
	R.set(ind(0)-1,old+1);
      }
      else {
	ivec ind = find(R == max(R));
	if (ind(0) == R.length()-1) {
	  R = concat(R,0);
	  Ri = linspace(0,length(R)-1,length(R));
	}
	int old = R(ind(0));
	R.set(ind(0),old-1);
	old = R(ind(0)+1);
	R.set(ind(0)+1,old+1);
      }
    }

    C = concat(C, zeros_i(Nmax-length(C)));
    R = concat(R, zeros_i(Nmax-length(R)));

    // -------------------

    if (method=="rand") {
      generate_random_H(C,R,options);
    } else { 
      it_error("not implemented"); 
    };
    
  }

  void LDPC_Parity_Matrix::generate_random_H(const ivec C, 
					     const ivec R, 
					     const ivec cycopt)
  {
    // Method based on random permutation. Attempts to avoid placing new
    // edges so that cycles are created. More aggressive cycle avoidance
    // for degree-2 nodes. EGL January 2007.

    initialize(sum(R),sum(C));
    // C, R: Target number of columns/rows with certain number of ones

    // compute number of edges
    int Ne=0;
    for (int i = 0;i < C.length();i++){
      for (int j = 0; j < C(i); j++) {
	for (int m=0; m<i; m++) Ne++;
      }
    }
 
    // compute connectivity matrix
    ivec vcon(Ne);
    ivec ccon(Ne); 
    ivec vd(nvar);
    ivec cd(ncheck); 
    int k=0;
    int l=0;
    for (int i = 0;i < C.length();i++){
      for (int j = 0; j < C(i); j++) {
	for (int m=0; m<i; m++) {
	  vcon(k)=l;
	  vd(l)=i;
	  k++;
	}
	l++;
      }
    }
    k=0;
    l=0;
    for (int i = 0;i < R.length();i++){
      for (int j = 0; j < R(i); j++) {
	for (int m=0; m<i; m++) {
	  ccon(k)=l;
	  cd(l)=i;
	  k++;
	} 
	l++;
      }
    }
    it_assert(k==Ne,"C/R mismatch");
  
    // compute random permutations
    ivec ind = sort_index(randu(Ne));   
    ivec cp = sort_index(randu(nvar));
    ivec rp = sort_index(randu(ncheck));

    // set the girth goal for various variable node degrees
    ivec Laim=zeros_i(Nmax);
    for (int i=0; i<length(cycopt); i++) {
      Laim(i+2)=cycopt(i);
    }
    for (int i=length(cycopt); i<Nmax-2; i++) {
      Laim(i+2)=cycopt(length(cycopt)-1);
    }
    cout << "Running with Laim=" << Laim.left(25) << endl;

    int failures=0;    
    const int Max_attempts=100;
    const int apcl=10;      // attempts before reducing girth target
    for (int k=0; k<Ne; k++) {
      const int el=Ne-k-2;
      if (k%250==0) { 
	cout << "Processing edge: " << k << " out of " << Ne 
	     << ". Variable node degree: " << vd(vcon(k)) 
	     << ". Girth target: " << Laim(vd(vcon(k))) 
	     << ". Accumulated failures: " << failures << endl ; 
      }
      const int c=cp(vcon(k));    
      int L= Laim(vd(vcon(k)));
      int attempt=0;
      while (true) {
	if (attempt>0 && attempt%apcl==0 && L>=6) { L-=2; };
	int r=rp(ccon(ind(k)));
	if (get(r,c)) { // double edge
	  // set(r,c,0);  
	  if (el>0) {
// 	    int t=k+1+rand()%el;
	    int t=k+1+randi(0,el-1);
	    int x=ind(t);  
	    ind(t)=ind(k);
	    ind(k)=x;  
	    attempt++;
	    if (attempt==Max_attempts) { 
	      failures++;
	      break; 
	    }
	  } else {  // almost at the last edge
	    break; 
	  }
	} else {
	  set(r,c,1);
	  if (L>0) { // attempt to avoid cycles
	    if (check_connectivity(r,c,r,c,0,L)>=0) {   // detected cycle
	      set(r,c,0);
	      if (el>0) {
		// make a swap in the index permutation
// 		int t=k+1+rand()%el;
		int t=k+1+randi(0,el-1);
		int x=ind(t);  
		ind(t)=ind(k);
		ind(k)=x;  
		attempt++;
		if (attempt==Max_attempts) {  // give up
		  failures++;
		  set(r,c,1);		
		  break;
		}
	      } else {  // no edges left
		set(r,c,1);		
		break; 
	      }
	    } else {
	      break;
	    }
	  } else {
	    break;
	  }
	}
      }
    }
  
    display_stats();
  }
  
  int LDPC_Parity_Matrix::check_connectivity(int from_i, int from_j, int to_i, 
				    int to_j, int godir, int L ) const
  {
    int i, j, result;
  
    if (L<0) {           // unable to reach coordinate with given L 
      return (-3); 
    }
    
    // check if reached destination 
    if ((from_i==to_i) && (from_j==to_j) && (godir!=0)) {  
      return L;
    }
  
    if (get(from_i,from_j)==0) {  // meaningless search 
      return (-2); 
    } 
  
    if (L==2) {      // Treat this case separately for efficiency
      if (godir==2) { // go horizontally
	if (get(from_i,to_j)==1) { return 0; }
      }
      if (godir==1) { // go vertically
	if (get(to_i,from_j)==1) { return 0; }
      }
      return (-3);
    }
  
    if ((godir==1) || (godir==0)) {   // go vertically
      ivec cj = index_nonzeros(get_col(from_j));
      for (i=0; i<length(cj); i++) {
	if (cj(i)!=from_i) {
	  result = check_connectivity(cj(i),from_j,to_i,to_j,2,L-1);
	  if (result>=0) {
	    return (result);
	  }
	}
      }
    }
  
    if (godir==2) {   // go horizontally
      ivec ri = index_nonzeros(get_row(from_i));
      for (j=0; j<length(ri); j++) {
	if (ri(j)!=from_j) {
	  result = check_connectivity(from_i,ri(j),to_i,to_j,1,L-1);
	  if (result>=0) {
	    return (result);
	  }
	}
      }
    }  
  
    return (-1);    
  };

  int LDPC_Parity_Matrix::check_for_cycles(int L) const
  {
    // looking for odd length cycles does not make sense
    if ((L&1)==1) { return (-1); } 
    if (L==0) { return (-4); } 

    int cycles=0;
    for (int i=0; i<nvar; i++) {
      ivec ri = index_nonzeros(get_col(i));
      for (int j=0; j<length(ri); j++) {
	if (check_connectivity(ri(j),i,ri(j),i,0,L)>=0)	{
	  cycles++;
	}
      }
    }
    return cycles;
  };

  // ---------------------------------------------------------------------------
  // LDPC_Generator_Matrix
  // ---------------------------------------------------------------------------

  LDPC_Generator_Matrix& LDPC_Generator_Matrix::operator=(const LDPC_Generator_Matrix &X)
  {
    G=X.G;
    type=X.type;
    nvar=X.nvar;
    ncheck=X.ncheck;
    return *this;
  };

  LDPC_Generator_Matrix::LDPC_Generator_Matrix(LDPC_Parity_Matrix &H, std::string t, 
					       ivec avoid_cols) 
  {
    it_assert(t=="systematic","LDPC_Generator_Matrix(): unsupported format");

    // create systematic generator matrix
    type=1;
    nvar=H.get_nvar();
    ncheck=H.get_ncheck();
    // systematic part is never saved explicitly

    // create dense representation of parity check matrix
    GF2mat Hd(H.H);  

    // take the columns in random order, but the ones to avoid at last
    vec col_importance = randu(nvar);
    for (int i=0; i<length(avoid_cols); i++) {
      col_importance(avoid_cols(i)) = (-col_importance(avoid_cols(i)));
    }
    ivec col_order = sort_index(-col_importance);

    GF2mat P1; //(ncheck,nvar-ncheck);      // non-invertible part
    GF2mat P2; //(ncheck,ncheck);           // invertible part
 
    cout << "Computing generator matrix..." << endl;

    int j1=0, j2=0;
    int rank;
    ivec perm;
    GF2mat T, U;
    for (int k=0; k<nvar; k++) {
      if (j1>=nvar-ncheck) {
	cout << "Unable to obtain enough independent columns." << endl;
	it_error("LDPC_Gmat(): internal error");
      }
      bvec c = Hd.get_col(col_order(k));
      if (j2==0)  {       // first column in P2 is number col_order(k)
	P2 = GF2mat(c);
	rank = P2.T_fact(T,U,perm);
	j2++;
      } else  {
	//      bvec c = Hd.get_col(col_order(k));
	if (j2<ncheck) {
	  bool success =  P2.T_fact_update_addcol(T,U,perm,c);
	  if (success) {
	    P2 = P2.concatenate_horizontal(c);
	    j2++;
	    continue;
	  } 
	} 
	if (j1==0) {
	  P1 = GF2mat(c); 
	} else {
	  P1 = P1.concatenate_horizontal(c);
	}
	j1++;
      }
    }

    cout << "Rank of parity check matrix: " << j2 << endl;

    // compute the systematic part of the generator matrix
    G = (P2.inverse()*P1).transpose();

    // permute the columns of the parity check matrix
    GF2mat P = P1.concatenate_horizontal(P2);
    H=LDPC_Parity_Matrix(ncheck,nvar);
    // brute force copy from P to H
    for (int i=0; i<ncheck; i++) {
      for (int j=0; j<nvar; j++) {
 	if (P.get(i,j)) {
 	  H.set(i,j,1);
 	}
      }
    }
  
    // check that the result was correct
    it_assert0((GF2mat(H.H)*(gf2dense_eye(nvar-ncheck)
			     .concatenate_horizontal(G)).transpose())
	       .is_zero(),
	       "LDPC_Gmat::LDPC_Gmat() result was incorrect");

    G=G.transpose();  // store the generator matrix in transposed form
    cout << "LDPC_Generator_Matrix() done." << endl;
  }


  // ---------------------------------------------------------------------------
  // LDPC_Code
  // ---------------------------------------------------------------------------

  LDPC_Code::LDPC_Code()
  { 
    H_is_defined=0;
    setup_decoder();
  } 

  LDPC_Code::LDPC_Code(const LDPC_Parity_Matrix &H)
  {
    set_code(H);
  }

  void LDPC_Code::set_code(const LDPC_Parity_Matrix &H)
  {
    H_is_defined=0;
    decoder_parameterization(H);
    G=LDPC_Generator_Matrix();  // clear the generator matrix

    setup_decoder();
  }

  void LDPC_Code::setup_decoder(std::string method_in, ivec options_in, 
				const LLR_calc_unit lcalc)
  {
    dec_method=method_in;
    dec_options=options_in;
    llrcalc=lcalc;   
    mcv.set_size(max(sumX2)*ncheck);
    mvc.set_size(nvar*max(sumX1));
  }

  LDPC_Code::LDPC_Code(const LDPC_Parity_Matrix &Hmat, 
		       const LDPC_Generator_Matrix &Gi)
  {
    set_code(Hmat,Gi);
  }
  
  void LDPC_Code::set_code(const LDPC_Parity_Matrix &Hmat, 
			   const LDPC_Generator_Matrix &Gi)
  {
    H_is_defined=0;
    decoder_parameterization(Hmat);
    G=Gi;
    
    // check that the generator and parity check matrices match
    GF2mat G1 = gf2dense_eye(nvar-ncheck)  
      .concatenate_horizontal(G.G.transpose());
    it_assert((GF2mat(Hmat.H)*G1.transpose()).is_zero(),
	      "LDPC_Code(): H and G matrices do not match");
    
//     mcv.set_size(max(sumX2)*ncheck);
//     mvc.set_size(nvar*max(sumX1));
    setup_decoder();
  }

  LDPC_Code::LDPC_Code(std::string filename)
  {
    set_code(filename);
  }

  void LDPC_Code::set_code(std::string filename)
  {
    cout << "LDPC_Code::LDPC_Code(): Loading LDPC codec from " 
	 << filename << endl;  
    it_ifile f(filename);
    int fileversion;
    f >> Name("Fileversion") >> fileversion;
    it_assert(fileversion==1,"unsupported format");
    f >> Name("H_is_defined") >> H_is_defined;
    f >> Name("nvar") >> nvar;
    f >> Name("ncheck") >> ncheck;
    f >> Name("C") >> C;
    f >> Name("V") >> V;
    f >> Name("sumX1") >> sumX1;
    f >> Name("sumX2") >> sumX2;
    f >> Name("iind") >> iind;
    f >> Name("jind") >> jind;
    f >> Name("G_type") >> G.type;
    switch (G.type) {
    case 0: { // no generator defined
      break;
    }
    case 1: { // systematic matrix
      f >> Name("G") >> G.G;
      // check integretity of generator matrix
      for (int i=0; i<G.G.cols(); i++) {
	bvec ci = G.G.get_col(i);
	ci = concat(zeros_b(G.G.cols()-1-i),ci);
	ci = concat(bin(1),ci);
	ci = concat(zeros_b(i),ci);
	it_assert(syndrome_check(ci),"LDPC_Code(): H and G do not match");
      }
      cout << "LDPC_Code::LDPC_Code(): Passed cross-check check." << endl;
      break;
    }
    default:
      it_error("LDPC_Code::LDPC_Code(): unsupported generator matrix format");
    }
    
    f.close();

    mcv.set_size(max(sumX2)*ncheck);
    mvc.set_size(nvar*max(sumX1));
    setup_decoder();

    cout << "LDPC_Code::LDPC_Code(): Successfully loaded LDPC codec from " 
	 << filename << endl;  
  }

  void LDPC_Code::save_to_file(std::string filename) const
  {
    cout << "Saving LDPC codec to " << filename << endl;
    it_file f;
    f.open(filename,true);
    f << Name("Fileversion") << 1;
    f << Name("H_is_defined") << H_is_defined;
    f << Name("nvar") << nvar;
    f << Name("ncheck") << ncheck;
    f << Name("C") << C;
    f << Name("V") << V;
    f << Name("sumX1") << sumX1;
    f << Name("sumX2") << sumX2;
    f << Name("iind") << iind;
    f << Name("jind") << jind;
    f << Name("G_type") << G.type;
    switch (G.type) {
    case 0: { // no generator
      break;
    }
    case 1: { // systematic
      f << Name("G") << G.G;
      break;
    }
    default:
      it_error("LDPC_Code::save_to_file(): unsupported G format");
    }

    f.close();

    cout << "LDPC_Code::save_to_file(): "
	 << " successfully saved LDPC codec to " << filename << endl;  
  }

  void LDPC_Code::decoder_parameterization(const LDPC_Parity_Matrix &Hmat)
  {
    // copy basic parameters
    sumX1 = Hmat.get_sumX1();
    sumX2 = Hmat.get_sumX2();
    nvar = Hmat.get_nvar();
    ncheck = Hmat.get_ncheck();
    int cmax = max(sumX1);
    int vmax = max(sumX2);

    // decoder parameterization
    V = zeros_i(ncheck*vmax);
    C = zeros_i(cmax*nvar);
    jind = zeros_i(ncheck*vmax);
    iind = zeros_i(nvar*cmax);
    
    cout << "Computing decoder parameterization. Phase 1." << endl;
    for (int i=0; i<nvar; i++) {
      ivec coli = Hmat.index_nonzeros(Hmat.get_col(i));
      for (int j0=0; j0<length(coli); j0++) {
	C(j0+cmax*i) = coli(j0);
      }
    }

    cout << "Computing decoder parameterization. Phase 2." << endl;
    for (int j=0; j<ncheck; j++) {
      ivec rowj = Hmat.index_nonzeros(Hmat.get_row(j));
      for (int i0=0; i0<length(rowj); i0++) {
	V(j+ncheck*i0) = rowj(i0);
      }
    }   

    cout << "Computing decoder parameterization. Phase 3." << endl;
    for (int j=0; j<ncheck; j++) {
      for (int ip=0; ip<sumX2(j); ip++) {
	int vip = V(j+ip*ncheck);
	int k=0;
	while (1==1)  { 
	  if (C(k+vip*cmax)==j) {
	    break;
	  }
	  k++; 
	}
	jind(j+ip*ncheck) = vip+k*nvar;
      }
    }

    cout << "Computing decoder parameterization. Phase 4." << endl;
    for (int i=0; i<nvar; i++) {
      for (int jp=0; jp<sumX1(i); jp++) {
	int cjp = C(jp+i*cmax);
	int k=0;
	while (1==1) {
	  if (V(cjp+k*ncheck)==i) {break; }
	  k++;
	}
	iind(i+jp*nvar) = cjp+k*ncheck;
      }
    }

    H_is_defined=1;
  }

  std::ostream &operator<<(std::ostream &os, LDPC_Code &C)
  {
    ivec rdeg = zeros_i(max(C.sumX2)+1);
    for (int i=0; i<C.ncheck; i++)     {
      rdeg(C.sumX2(i))++;
    }  

    ivec cdeg = zeros_i(max(C.sumX1)+1);
    for (int j=0; j<C.nvar; j++)     {
      cdeg(C.sumX1(j))++;
    }
    
    os << "---------- LDPC codec -----------------" << std::endl;
    os << "Nvar=" << C.get_nvar() << "  Ncheck=" << C.get_ncheck() 
       << "  rate=" << C.get_rate() << std::endl;
    os << "Column degrees (node perspective): " << cdeg << std::endl;
    os << "Row degrees (node perspective): " << rdeg << std::endl;
    os << "---------------------------------------" << std::endl;
    os << "Decoder parameters: dec_method=" << C.dec_method 
       << ". Decoder options: " << C.dec_options << std::endl;
    os << C.llrcalc << std::endl;
    return os;
  }

  bool LDPC_Code::syndrome_check(const bvec &x) const
  {
    QLLRvec llr=1-2*to_ivec(x);
    return syndrome_check(llr);
  }
  
  bool LDPC_Code::syndrome_check(const QLLRvec &LLR) const    
  {
    // Please note the IT++ convention that a sure zero corresponds to
    // LLR=+infinity
    int i,j,synd,vi;
  
    for (j=0; j<ncheck; j++) {
      synd = 0;
      int vind = j; // tracks j+i*ncheck
      for (i=0; i<sumX2(j); i++) {
	vi = V(vind);
	if (LLR(vi)<0) {
	  synd++;
	}
	vind += ncheck;
      }
      if ((synd&1)==1) { 
	return false;  // codeword is invalid
      }    
    }          
    return true;   // codeword is valid 
  };


  int LDPC_Code::bp_decode(const QLLRvec &LLRin,   QLLRvec &LLRout) 
  {
    // Note the IT++ convention that a sure zero corresponds to
    // LLR=+infinity
    it_assert(H_is_defined,"there is no parity check matrix");
    it_assert1(length(LLRin)==nvar 
	       && length(sumX1)==nvar 
	       && length(sumX2)==ncheck,
	       "LDPC_Code::bp_decode() wrong input dimensions");
    LLRout.set_size(length(LLRin));

    const int niter=dec_options(0);
    const int psc=dec_options(1);
    const int pisc=dec_options(2);
  
    if (pisc==1) {
      if (syndrome_check(LLRin)) {
	LLRout = LLRin;
	return 0; 
      }
    }

    // initial step
    for (int i=0; i<nvar; i++) {
      int index = i;
      for (int j=0; j<sumX1(i); j++) {
	mvc[index] = LLRin(i); 
	index += nvar;
      }
    }

    bool is_valid_codeword=false;
    int iter =0;
    do {
      iter++;  
      if (nvar>=100000) { cout << "."; cout.flush(); }
      // --------- Step 1: check to variable nodes ----------
      for (int j=0; j<ncheck; j++) {
	switch (sumX2(j)) {
	case 0: it_error("LDPC_Code::bp_decode(): sumX2(j)=0");
	case 1: it_error("LDPC_Code::bp_decode(): sumX2(j)=1");
	case 2: {
	  mcv[j+ncheck]=mvc[jind[j]];
	  mcv[j]=mvc[jind[j+ncheck]];
	  break;
	}
	case 3: {
	  int j0=j;
	  QLLR m0=mvc[jind[j0]];
	  int j1=j0+ncheck;
	  QLLR m1=mvc[jind[j1]];
	  int j2=j1+ncheck;
	  QLLR m2=mvc[jind[j2]];
	  mcv[j0]=llrcalc.Boxplus(m1,m2);
	  mcv[j1]=llrcalc.Boxplus(m0,m2);
	  mcv[j2]=llrcalc.Boxplus(m0,m1);
	  break;
	}
	case 4: {
	  int j0=j;
	  QLLR m0=mvc[jind[j0]];
	  int j1=j0+ncheck;
	  QLLR m1=mvc[jind[j1]];
	  int j2=j1+ncheck;
	  QLLR m2=mvc[jind[j2]];
	  int j3=j2+ncheck;
	  QLLR m3=mvc[jind[j3]];
	  QLLR m01=llrcalc.Boxplus(m0,m1);
	  QLLR m23=llrcalc.Boxplus(m2,m3);
	  mcv[j0]=llrcalc.Boxplus(m1,m23);
	  mcv[j1]=llrcalc.Boxplus(m0,m23);
	  mcv[j2]=llrcalc.Boxplus(m01,m3);
	  mcv[j3]=llrcalc.Boxplus(m01,m2);
	  break;
	} 
	case 5: {
	  int j0=j;
	  QLLR m0=mvc[jind[j0]];
	  int j1=j0+ncheck;
	  QLLR m1=mvc[jind[j1]];
	  int j2=j1+ncheck;
	  QLLR m2=mvc[jind[j2]];
	  int j3=j2+ncheck;
	  QLLR m3=mvc[jind[j3]];
	  int j4=j3+ncheck;
	  QLLR m4=mvc[jind[j4]];
	  QLLR m01=llrcalc.Boxplus(m0,m1);
	  QLLR m02=llrcalc.Boxplus(m01,m2);
	  QLLR m34=llrcalc.Boxplus(m3,m4);
	  QLLR m24=llrcalc.Boxplus(m2,m34);
	  mcv[j0]=llrcalc.Boxplus(m1,m24);
	  mcv[j1]=llrcalc.Boxplus(m0,m24);
	  mcv[j2]=llrcalc.Boxplus(m01,m34);
	  mcv[j3]=llrcalc.Boxplus(m02,m4);
	  mcv[j4]=llrcalc.Boxplus(m02,m3);
	  break;
	}
	case 6: {
	  int j0=j;
	  QLLR m0=mvc[jind[j0]];
	  int j1=j0+ncheck;
	  QLLR m1=mvc[jind[j1]];
	  int j2=j1+ncheck;
	  QLLR m2=mvc[jind[j2]];
	  int j3=j2+ncheck;
	  QLLR m3=mvc[jind[j3]];
	  int j4=j3+ncheck;
	  QLLR m4=mvc[jind[j4]];
	  int j5=j4+ncheck;
	  QLLR m5=mvc[jind[j5]];
	  QLLR m01=llrcalc.Boxplus(m0,m1);
	  QLLR m23=llrcalc.Boxplus(m2,m3);
	  QLLR m45=llrcalc.Boxplus(m4,m5);
	  QLLR m03=llrcalc.Boxplus(m01,m23);
	  QLLR m25=llrcalc.Boxplus(m23,m45);
	  QLLR m0145=llrcalc.Boxplus(m01,m45);
	  mcv[j0]=llrcalc.Boxplus(m1,m25);
	  mcv[j1]=llrcalc.Boxplus(m0,m25);
	  mcv[j2]=llrcalc.Boxplus(m0145,m3);
	  mcv[j3]=llrcalc.Boxplus(m0145,m2);
	  mcv[j4]=llrcalc.Boxplus(m03,m5);
	  mcv[j5]=llrcalc.Boxplus(m03,m4);
	  break;
	}	
	case 7: {
	  int j0=j;
	  QLLR m0=mvc[jind[j0]];
	  int j1=j0+ncheck;
	  QLLR m1=mvc[jind[j1]];
	  int j2=j1+ncheck;
	  QLLR m2=mvc[jind[j2]];
	  int j3=j2+ncheck;
	  QLLR m3=mvc[jind[j3]];
	  int j4=j3+ncheck;
	  QLLR m4=mvc[jind[j4]];
	  int j5=j4+ncheck;
	  QLLR m5=mvc[jind[j5]];
	  int j6=j5+ncheck;
	  QLLR m6=mvc[jind[j6]];
	  QLLR m01=llrcalc.Boxplus(m0,m1);
	  QLLR m23=llrcalc.Boxplus(m2,m3);
	  QLLR m45=llrcalc.Boxplus(m4,m5);
	  QLLR m46=llrcalc.Boxplus(m45,m6);
	  QLLR m03=llrcalc.Boxplus(m01,m23);
	  QLLR m25=llrcalc.Boxplus(m23,m45);
	  QLLR m26=llrcalc.Boxplus(m25,m6);
	  QLLR m04=llrcalc.Boxplus(m03,m4);
	  mcv[j0]=llrcalc.Boxplus(m26,m1);
	  mcv[j1]=llrcalc.Boxplus(m26,m0);
	  mcv[j2]=llrcalc.Boxplus(m01,llrcalc.Boxplus(m3,m46));
	  mcv[j3]=llrcalc.Boxplus(m2,llrcalc.Boxplus(m01,m46));
	  mcv[j4]=llrcalc.Boxplus(m6,llrcalc.Boxplus(m03,m5));
	  mcv[j5]=llrcalc.Boxplus(m6,m04);
	  mcv[j6]=llrcalc.Boxplus(m5,m04);
	  break;
	}
	case 8: {
	  int j0=j;
	  QLLR m0=mvc[jind[j0]];
	  int j1=j0+ncheck;
	  QLLR m1=mvc[jind[j1]];
	  int j2=j1+ncheck;
	  QLLR m2=mvc[jind[j2]];
	  int j3=j2+ncheck;
	  QLLR m3=mvc[jind[j3]];
	  int j4=j3+ncheck;
	  QLLR m4=mvc[jind[j4]];
	  int j5=j4+ncheck;
	  QLLR m5=mvc[jind[j5]];
	  int j6=j5+ncheck;
	  QLLR m6=mvc[jind[j6]];
	  int j7=j6+ncheck;
	  QLLR m7=mvc[jind[j7]];
	  QLLR m01=llrcalc.Boxplus(m0,m1);
	  QLLR m23=llrcalc.Boxplus(m2,m3);
	  QLLR m45=llrcalc.Boxplus(m4,m5);
	  QLLR m67=llrcalc.Boxplus(m6,m7);
	  QLLR m47=llrcalc.Boxplus(m45,m67);
	  QLLR m03=llrcalc.Boxplus(m01,m23);
	  QLLR m25=llrcalc.Boxplus(m23,m45);
	  mcv[j0]=llrcalc.Boxplus(m67,llrcalc.Boxplus(m1,m25));
	  mcv[j1]=llrcalc.Boxplus(m67,llrcalc.Boxplus(m0,m25));
	  mcv[j2]=llrcalc.Boxplus(m3,llrcalc.Boxplus(m01,m47));
	  mcv[j3]=llrcalc.Boxplus(m2,llrcalc.Boxplus(m01,m47));
	  mcv[j4]=llrcalc.Boxplus(m67,llrcalc.Boxplus(m03,m5));
	  mcv[j5]=llrcalc.Boxplus(m67,llrcalc.Boxplus(m03,m4));
	  mcv[j6]=llrcalc.Boxplus(m45,llrcalc.Boxplus(m03,m7));
	  mcv[j7]=llrcalc.Boxplus(m03,llrcalc.Boxplus(m45,m6));
	  break;
	}	
	case 9: {
	  int j0=j;
	  QLLR m0=mvc[jind[j0]];
	  int j1=j0+ncheck;
	  QLLR m1=mvc[jind[j1]];
	  int j2=j1+ncheck;
	  QLLR m2=mvc[jind[j2]];
	  int j3=j2+ncheck;
	  QLLR m3=mvc[jind[j3]];
	  int j4=j3+ncheck;
	  QLLR m4=mvc[jind[j4]];
	  int j5=j4+ncheck;
	  QLLR m5=mvc[jind[j5]];
	  int j6=j5+ncheck;
	  QLLR m6=mvc[jind[j6]];
	  int j7=j6+ncheck;
	  QLLR m7=mvc[jind[j7]];
	  int j8=j7+ncheck;
	  QLLR m8=mvc[jind[j8]];
	  QLLR m01=llrcalc.Boxplus(m0,m1);
	  QLLR m23=llrcalc.Boxplus(m2,m3);
	  QLLR m45=llrcalc.Boxplus(m4,m5);
	  QLLR m67=llrcalc.Boxplus(m6,m7);
	  QLLR m68=llrcalc.Boxplus(m67,m8);
	  QLLR m03=llrcalc.Boxplus(m01,m23);
	  QLLR m25=llrcalc.Boxplus(m23,m45);
	  QLLR m05=llrcalc.Boxplus(m03,m45);
	  mcv[j0]=llrcalc.Boxplus(m68,llrcalc.Boxplus(m1,m25));
	  mcv[j1]=llrcalc.Boxplus(m68,llrcalc.Boxplus(m0,m25));
	  mcv[j2]=llrcalc.Boxplus(llrcalc.Boxplus(m01,m68),
				  llrcalc.Boxplus(m3,m45));
	  mcv[j3]=llrcalc.Boxplus(llrcalc.Boxplus(m01,m68),
				  llrcalc.Boxplus(m2,m45));
	  mcv[j4]=llrcalc.Boxplus(m68,llrcalc.Boxplus(m03,m5));
	  mcv[j5]=llrcalc.Boxplus(m68,llrcalc.Boxplus(m03,m4));
	  mcv[j6]=llrcalc.Boxplus(llrcalc.Boxplus(m7,m8),m05);
	  mcv[j7]=llrcalc.Boxplus(llrcalc.Boxplus(m05,m6),m8);
	  mcv[j8]=llrcalc.Boxplus(m05,m67);
	  break;
	}	
	case 10: {
	  int j0=j;
	  QLLR m0=mvc[jind[j0]];
	  int j1=j0+ncheck;
	  QLLR m1=mvc[jind[j1]];
	  int j2=j1+ncheck;
	  QLLR m2=mvc[jind[j2]];
	  int j3=j2+ncheck;
	  QLLR m3=mvc[jind[j3]];
	  int j4=j3+ncheck;
	  QLLR m4=mvc[jind[j4]];
	  int j5=j4+ncheck;
	  QLLR m5=mvc[jind[j5]];
	  int j6=j5+ncheck;
	  QLLR m6=mvc[jind[j6]];
	  int j7=j6+ncheck;
	  QLLR m7=mvc[jind[j7]];
	  int j8=j7+ncheck;
	  QLLR m8=mvc[jind[j8]];
	  int j9=j8+ncheck;
	  QLLR m9=mvc[jind[j9]];
	  QLLR m01=llrcalc.Boxplus(m0,m1);
	  QLLR m23=llrcalc.Boxplus(m2,m3);
	  QLLR m03=llrcalc.Boxplus(m01,m23);
	  QLLR m45=llrcalc.Boxplus(m4,m5);
	  QLLR m67=llrcalc.Boxplus(m6,m7);
	  QLLR m89=llrcalc.Boxplus(m8,m9);
	  QLLR m69=llrcalc.Boxplus(m67,m89);
	  QLLR m25=llrcalc.Boxplus(m23,m45);
	  QLLR m05=llrcalc.Boxplus(m03,m45);
	  QLLR m07=llrcalc.Boxplus(m05,m67);
	  mcv[j0]=llrcalc.Boxplus(m69,llrcalc.Boxplus(m1,m25));
	  mcv[j1]=llrcalc.Boxplus(m69,llrcalc.Boxplus(m0,m25));
	  mcv[j2]=llrcalc.Boxplus(llrcalc.Boxplus(m01,m69),
				  llrcalc.Boxplus(m3,m45));
	  mcv[j3]=llrcalc.Boxplus(llrcalc.Boxplus(m01,m69),
				  llrcalc.Boxplus(m2,m45));
	  mcv[j4]=llrcalc.Boxplus(m69,llrcalc.Boxplus(m03,m5));
	  mcv[j5]=llrcalc.Boxplus(m69,llrcalc.Boxplus(m03,m4));
	  mcv[j6]=llrcalc.Boxplus(llrcalc.Boxplus(m7,m89),m05);
	  mcv[j7]=llrcalc.Boxplus(llrcalc.Boxplus(m05,m6),m89);
	  mcv[j8]=llrcalc.Boxplus(m07,m9);
	  mcv[j9]=llrcalc.Boxplus(m07,m8);
	  break;
	}	
	case 11: {
	  int j0=j;
	  QLLR m0=mvc[jind[j0]];
	  int j1=j0+ncheck;
	  QLLR m1=mvc[jind[j1]];
	  int j2=j1+ncheck;
	  QLLR m2=mvc[jind[j2]];
	  int j3=j2+ncheck;
	  QLLR m3=mvc[jind[j3]];
	  int j4=j3+ncheck;
	  QLLR m4=mvc[jind[j4]];
	  int j5=j4+ncheck;
	  QLLR m5=mvc[jind[j5]];
	  int j6=j5+ncheck;
	  QLLR m6=mvc[jind[j6]];
	  int j7=j6+ncheck;
	  QLLR m7=mvc[jind[j7]];
	  int j8=j7+ncheck;
	  QLLR m8=mvc[jind[j8]];
	  int j9=j8+ncheck;
	  QLLR m9=mvc[jind[j9]];
	  int j10=j9+ncheck;
	  QLLR m10=mvc[jind[j10]];
	  QLLR m01=llrcalc.Boxplus(m0,m1);
	  QLLR m23=llrcalc.Boxplus(m2,m3);
	  QLLR m03=llrcalc.Boxplus(m01,m23);
	  QLLR m45=llrcalc.Boxplus(m4,m5);
	  QLLR m67=llrcalc.Boxplus(m6,m7);
	  QLLR m89=llrcalc.Boxplus(m8,m9);
	  QLLR m69=llrcalc.Boxplus(m67,m89);
	  QLLR m6_10=llrcalc.Boxplus(m69,m10);
	  QLLR m25=llrcalc.Boxplus(m23,m45);
	  QLLR m05=llrcalc.Boxplus(m03,m45);
	  QLLR m07=llrcalc.Boxplus(m05,m67);
	  QLLR m8_10=llrcalc.Boxplus(m89,m10);
	  mcv[j0]=llrcalc.Boxplus(m6_10,llrcalc.Boxplus(m1,m25));
	  mcv[j1]=llrcalc.Boxplus(m6_10,llrcalc.Boxplus(m0,m25));
	  mcv[j2]=llrcalc.Boxplus(llrcalc.Boxplus(m01,m6_10),
				  llrcalc.Boxplus(m3,m45));
	  mcv[j3]=llrcalc.Boxplus(llrcalc.Boxplus(m01,m6_10),
				  llrcalc.Boxplus(m2,m45));
	  mcv[j4]=llrcalc.Boxplus(m6_10,llrcalc.Boxplus(m03,m5));
	  mcv[j5]=llrcalc.Boxplus(m6_10,llrcalc.Boxplus(m03,m4));
	  mcv[j6]=llrcalc.Boxplus(llrcalc.Boxplus(m7,m8_10),m05);
	  mcv[j7]=llrcalc.Boxplus(llrcalc.Boxplus(m05,m6),m8_10);
	  mcv[j8]=llrcalc.Boxplus(m10,llrcalc.Boxplus(m07,m9));
	  mcv[j9]=llrcalc.Boxplus(m10,llrcalc.Boxplus(m07,m8));
	  mcv[j10]=llrcalc.Boxplus(m07,m89);
	  break;
	}	
	case 12: {
	  int j0=j;
	  QLLR m0=mvc[jind[j0]];
	  int j1=j0+ncheck;
	  QLLR m1=mvc[jind[j1]];
	  int j2=j1+ncheck;
	  QLLR m2=mvc[jind[j2]];
	  int j3=j2+ncheck;
	  QLLR m3=mvc[jind[j3]];
	  int j4=j3+ncheck;
	  QLLR m4=mvc[jind[j4]];
	  int j5=j4+ncheck;
	  QLLR m5=mvc[jind[j5]];
	  int j6=j5+ncheck;
	  QLLR m6=mvc[jind[j6]];
	  int j7=j6+ncheck;
	  QLLR m7=mvc[jind[j7]];
	  int j8=j7+ncheck;
	  QLLR m8=mvc[jind[j8]];
	  int j9=j8+ncheck;
	  QLLR m9=mvc[jind[j9]];
	  int j10=j9+ncheck;  
	  QLLR m10=mvc[jind[j10]];
	  int j11=j10+ncheck;
	  QLLR m11=mvc[jind[j11]];
	  QLLR m01=llrcalc.Boxplus(m0,m1);
	  QLLR m23=llrcalc.Boxplus(m2,m3);
	  QLLR m03=llrcalc.Boxplus(m01,m23);
	  QLLR m45=llrcalc.Boxplus(m4,m5);
	  QLLR m67=llrcalc.Boxplus(m6,m7);
	  QLLR m89=llrcalc.Boxplus(m8,m9);
	  QLLR m69=llrcalc.Boxplus(m67,m89);
	  QLLR m10_11=llrcalc.Boxplus(m10,m11);
	  QLLR m6_11=llrcalc.Boxplus(m69,m10_11);
	  QLLR m25=llrcalc.Boxplus(m23,m45);
	  QLLR m05=llrcalc.Boxplus(m03,m45);
	  QLLR m07=llrcalc.Boxplus(m05,m67);
	  QLLR m8_10=llrcalc.Boxplus(m89,m10);
	  mcv[j0]=llrcalc.Boxplus(m6_11,llrcalc.Boxplus(m1,m25));
	  mcv[j1]=llrcalc.Boxplus(m6_11,llrcalc.Boxplus(m0,m25));
	  mcv[j2]=llrcalc.Boxplus(llrcalc.Boxplus(m01,m6_11),
				  llrcalc.Boxplus(m3,m45));
	  mcv[j3]=llrcalc.Boxplus(llrcalc.Boxplus(m01,m6_11),
				  llrcalc.Boxplus(m2,m45));
	  mcv[j4]=llrcalc.Boxplus(m6_11,llrcalc.Boxplus(m03,m5));
	  mcv[j5]=llrcalc.Boxplus(m6_11,llrcalc.Boxplus(m03,m4));
	  mcv[j6]=llrcalc.Boxplus(m11,llrcalc.Boxplus(llrcalc.Boxplus(m7,m8_10),m05));
	  mcv[j7]=llrcalc.Boxplus(m11,llrcalc.Boxplus(llrcalc.Boxplus(m05,m6),m8_10));
	  mcv[j8]=llrcalc.Boxplus(m10_11,llrcalc.Boxplus(m07,m9));
	  mcv[j9]=llrcalc.Boxplus(m10_11,llrcalc.Boxplus(m07,m8));
	  mcv[j10]=llrcalc.Boxplus(llrcalc.Boxplus(m07,m89),m11);
	  mcv[j11]=llrcalc.Boxplus(llrcalc.Boxplus(m07,m89),m10);
	  break;
	}	
	default:  
	  it_error("check node degrees >12 not supported in this version");
	}  // switch statement
      }
    
      // step 2: variable to check nodes 
      for (int i=0; i<nvar; i++) {
	switch (sumX1(i)) {
	case 0: it_error("LDPC_Code::bp_decode(): sumX1(i)=0");
	case 1: {
	  // Degenerate case-should not occur for good coded. A lonely
	  // variable node only sends its incoming message
	  mvc[i] = LLRin(i); 
	  LLRout(i)=LLRin(i);
	  break;
	}
	case 2: {
	  QLLR m0=mcv[iind[i]];
	  int i1=i+nvar;
	  QLLR m1=mcv[iind[i1]];
	  mvc[i] = LLRin(i) + m1;
	  mvc[i1] = LLRin(i) + m0;
	  LLRout(i) = mvc[i1]+m1;
	  break;
	}
	case 3: {
	  int i0=i;
	  QLLR m0 = mcv[iind[i0]];
	  int i1 = i0+nvar;
	  QLLR m1 = mcv[iind[i1]];
	  int i2 = i1+nvar;
	  QLLR m2 = mcv[iind[i2]];
	  LLRout(i) = LLRin(i)+m0+m1+m2;
	  mvc[i0]=LLRout(i)-m0;
	  mvc[i1]=LLRout(i)-m1;
	  mvc[i2]=LLRout(i)-m2;
	  break;
	}
	case 4: {
	  int i0=i;
	  QLLR m0 = mcv[iind[i0]];
	  int i1 = i0+nvar;
	  QLLR m1 = mcv[iind[i1]];
	  int i2 = i1+nvar;
	  QLLR m2 = mcv[iind[i2]];
	  int i3 = i2+nvar;
	  QLLR m3 = mcv[iind[i3]];
	  LLRout(i)= LLRin(i)+m0+m1+m2+m3;
	  mvc[i0]=LLRout(i)-m0;
	  mvc[i1]=LLRout(i)-m1;
	  mvc[i2]=LLRout(i)-m2;
	  mvc[i3]=LLRout(i)-m3;
	  break;
	} 
	default:  	{ // differential update 
	  QLLR mvc_temp = LLRin(i);
	  int index_iind = i; // tracks i+jp*nvar
	  for (int jp=0; jp<sumX1(i); jp++) {
	    mvc_temp +=  mcv[iind[index_iind]];
	    index_iind += nvar;
	  }
	  LLRout(i) = mvc_temp;
	  index_iind = i;  // tracks i+j*nvar
	  for (int j=0; j<sumX1[i]; j++) {
	    mvc[index_iind] = mvc_temp - mcv[iind[index_iind]];
	    index_iind += nvar;
	  }
	}
	}
      }
      
      if (psc==1) {
	if (syndrome_check(LLRout)==1) {
	  is_valid_codeword=true;
	  break;
	}
      }
    } while  (iter<=niter);
   
    if (nvar>=100000) { cout << endl; }
    return (is_valid_codeword ? iter : -iter);
  }

  void LDPC_Code::decode(const vec &llrin, bvec &bitsout)
  {
    it_assert(dec_method=="bp","not implemented");
    it_assert(H_is_defined,"there is no parity check matrix");

    // decode with belief propagation
    QLLRvec qllrin=llrcalc.to_qllr(llrin);
    QLLRvec qllrout;
    bp_decode(qllrin,qllrout);
    bitsout = qllrout<0;
  }
  
  bvec LDPC_Code::decode(const vec &x)
  {
    it_assert(dec_method=="bp","not implemented");
    it_assert(H_is_defined,"there is no parity check matrix");

    // decode with belief propagation
    bvec b;
    decode(x,b);
    return b;
  }

  void LDPC_Code::encode(const bvec &bitsin, bvec &bitsout)
  {
    it_assert(G.type!=0,
	      "encode(): No G matrix, only zero codeword possible.");

    bitsout = zeros_b(nvar);
    switch (G.type) {
    case 1: {       // generator matrix is systematic
      for (int i=0; i<nvar-ncheck; i++) {
	bitsout(i)=bitsin(i);
      }
      bvec btemp=G.G*bitsin;  
      for (int i=0; i<ncheck; i++) {
	bitsout(i+nvar-ncheck) = btemp(i);
      }
      break;
    }
    default:
      it_error("LDPC_Code::encode(): unsupported G matrix type");
    }
    
    it_assert0(syndrome_check(bitsout),"LDPC_Code::encode()");
  } 

  bvec LDPC_Code::encode(const bvec &input)
  { 
    bvec result;
    encode(input,result);
    return result;
  }
}
