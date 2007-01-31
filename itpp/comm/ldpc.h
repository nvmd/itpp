/*!
 * \file 
 * \brief Implementation of a Low-Density Parity Check (LDPC) codec.
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

#ifndef LDPC_H
#define LDPC_H

#include <iostream>
#include <itpp/base/gf2mat.h>
#include <itpp/base/random.h>
#include <itpp/base/sort.h>
#include <itpp/comm/llr.h>
#include <itpp/comm/channel_code.h>

namespace itpp {
  
  // ---------------------------------------------------------------------------
  // LDPC_Parity_Matrix
  // ---------------------------------------------------------------------------
  
  /*!
    \brief LDPC Parity check matrix class

    This class provides functionality for design of and operations on
    parity check matrices for LDPC codes. Please refer to the tutorial
    \ref ldpc_gen_codes for extensive examples of code generation.
    
    Parity check matrices can be saved to a file by converting them to
    \c GF2mat_sparse_alist format. However, typically one will want to
    create a \c LDPC_Code from the matrix and save the codec instead.

    This class stores a parity check matrix as a sparse matrix.  The
    transpose of the matrix is also stored to enable efficient access
    to its rows.

    \author Erik G. Larsson. Contributions from Mattias Andersson.
  */
  class LDPC_Parity_Matrix
    {
      friend class LDPC_Generator_Matrix;
      friend class LDPC_Code;
      
    public:     
      // ------------ Constructors ------------------------------
      
      /*! \brief Default constructor

      Gives an empty parity check matrix of size 1*1 */
      LDPC_Parity_Matrix();

      //! Constructor, gives empty matrix of size ncheck*nvar
      LDPC_Parity_Matrix(int ncheck, int nvar);

      /*! \brief Obtain LDPC parity check matrix from a file

      \param filename Filename 
      \param format format of the
      matrix. Currently, only "alist" is supported (see 
      \c GF2mat_sparse_alist for definition).
      */
      LDPC_Parity_Matrix(std::string filename, std::string format);

      //! Constructor, from a \c GF2mat_sparse_alist object
      LDPC_Parity_Matrix(const GF2mat_sparse_alist &alist);

      // ----------- Basic access functions ---------------------

      //! Get the parity check matrix 
      GF2mat_sparse get_H() const { return H; } 

      //! Get a specific column from the matrix
      Sparse_Vec<bin> get_col(int r) const {   return H.get_col(r); }

      //! Get a specific row from the matrix
      Sparse_Vec<bin> get_row(int r) const {   return Ht.get_col(r); }

      //! Get the number of variable nodes (number of columns) 
      int get_nvar() const {
	it_assert0(H.cols()==nvar,"LDPC_Parity_Matrix internal error");
	it_assert0(Ht.rows()==nvar,"LDPC_Parity_Matrix internal error");    
	return nvar; 
      }

      //! Get the number of check nodes (number of rows)
      int get_ncheck() const {
	it_assert0(H.rows()==ncheck,"LDPC_Parity_Matrix internal error");
	it_assert0(Ht.cols()==ncheck,"LDPC_Parity_Matrix internal error");
	return ncheck;
      }

      //! Set element (i,j) of the parity check matrix to value
      void set(int i, int j, bin value);

      //! Get element (i,j) of the parity check matrix
      inline bin get(int i, int j) const {   
	it_assert0(H(i,j)==Ht(j,i),"LDPC_Parity_Matrix internal error");
	return H(i,j);  
      }

      //! Display some information about the matrix
      void display_stats() const;

      //! Get the coderate
      inline double get_rate() const { 
	return   (1.0 - static_cast<double>(ncheck) / nvar); 
      }
      
      //! Export matrix to \c GF2mat_sparse_alist format
      GF2mat_sparse_alist export_alist() const;
  
      // ------------ Parity matrix generation ---------------------

      /*! \brief Generate a (k,l) regular LDPC code. 

      This function prepares for the generation of a regular LDPC
      code.

      \param Nvar  number of variable nodes
      \param k number of ones per column
      \param l number of ones per row
      \param method See documentation of generate_irregular_ldpc()
      \param options See documentation of generate_irregular_ldpc()
      */
      void generate_regular_ldpc(int Nvar, int k, int l, 
				 std::string method="rand",
				 ivec options="200 6");

      /*! \brief Generate an irregular LDPC code.

      This function  generates an irregular LDPC code. 

      \param Nvar  number of variable nodes

      \param var_deg vector of variable degree distributions, from an
      edge perspective

      \param chk_deg vector of check degree distributions, from an
      edge perspective

      \param method currently the only supported method is "rand" 
      (see below)

      \param options Determines the level of matrix optimization.  
      
      The "rand" method generates a fully unstructured random
      matrix. Some stochastic optimization is performed to avoid
      cycles. This optimization is controlled via the parameter \c
      options.  In particular, the girth for the
      variable-node-degree-2 part of the graph can be controlled via
      the parameter options(0). The girth for the rest of the nodes
      can be controlled via the parameter options(1).  The recommended
      value for options is "200 6" for graphs of small size, and "100
      4" for large graphs. The value "0 0" means no
      optimization. Double edges are always avoided.
      */
      void generate_irregular_ldpc(int Nvar, vec var_deg, vec chk_deg, 
				   std::string method="rand",
				   ivec options="200 6");

      // ----------- Graph optimization ------------------------

      /*!   \brief Remove cycles (loops) from parity check matrix.

      This function implements the cycle removal algorithm presented
      by McGowan and Williamson at the IT workshop 2003. The maximum
      girth of the graph that will be attempted is L. The algorithm is
      bound to remove all loops of length L, insofar this is
      possible. I.e., it does not terminate until it is impossible to
      remove more cycles by swapping two edges.

      \note This algorithm can take a long time to run for large L
      or large graphs.

      The function returns the girth of the graph, i.e. the length of
      the shortest cycle. For example, a return value of 6 means that
      there are no 4-cycles.

      \param L target girth. For example, L=6 attempts to removes all
      4-cycles.
      */
      int cycle_removal_MGW(int L);
 
    protected:
      //! Initialize empty matrix of size ncheck*nvar
      void initialize(int ncheck, int nvar);

      //! Get sum over all rows (sum(X,1))
      ivec get_sumX1() const {  return sumX1;   }

      //! Get sum over all columns (sum(X,2))
      ivec get_sumX2() const {  return sumX2;   }  

      /*! \brief Check for cycles of length L 

      This function implements a recursive routine to find loops. The
      function is mainly a tool for testing and debugging more
      sophisticated functions for graph manipulation. The function
      returns the number of cycles found of length L or
      shorter. Cycles may be counted multiple times.

      \param L length of cycles to look for

      \note This function can be very slow for large matrices. It is
      mainly intended as a debugging aid.
      */
      int check_for_cycles(int L) const;

      /*! \brief Check for connectivity between nodes 
	
      This function examines whether the point (to_m, to_n) in the
      matrix can be reached from the point (from_m, from_n) using at
      most L steps. A recursive search is used. 

      The function can be used to search for cycles in the matrix. To
      search for a cycle of length L, set from_m=to_m and from_n=to_n,
      and godir=0.

      \param from_m starting coordinate, row number
      \param to_m goal coordinte, row number
      \param from_n starting coordinate, column number
      \param to_n goal coordinate, row number
      \param g direction: 1=start going vertically, 2=start 
      going horizontally
      \param L number of permitted steps
     
      Possible return values:

      -1 or -3 : destination unreachable, 

      -2       : meaningless search (started in a "0" point), 

      -4       : meaningless search

      >=0      : destination reached with certain number of steps left

      \note This function can be very slow depending on the nature of
      the matrix.

      Note that smaller cycles may appear as longer cycles when using
      this method. More specifically, suppose the method is run with a
      given L and there are cycles in the neighborhood of
      (from_m,from_n) of length L-2 or less, but which do not contain
      (from_m,from_n). These shorter cycles may then also be reported
      as a cycle of length L.  For example, if one of the immediate
      neighbors of (from_m,from_n) is part of a cycle of length 4 this
      method will report that (from_m,from_n) is part of a cycle of
      length 6, if run with L=6. However, if it is known that there
      are no cycles of length L-2 or smaller, and
      check_connectivity(from_m,from_n,from_m,from_n,0,L) returns a
      non-negative value, then one will know with certainty that the
      point (from_m,from_n) is part of a cycle of length L. (This
      behavior is inherent to the simple recursive search used.)
      */
      int check_connectivity(int from_m,int from_n,int to_m,int to_n,int g,int L) const;

    private:
      static const int Nmax=200;   // Maximum node degree class can handle

      GF2mat_sparse H;    // The parity check matrix 
      GF2mat_sparse Ht;   // The transposed parity check matrix
      int nvar, ncheck;   // Matrix dimensions
      ivec sumX1, sumX2;  // Actual number of ones in each column and row
     
      // Generate random parity matrix
      void generate_random_H(const ivec C, const ivec R, const ivec cycopt);
      
      /* This help function returns the nonzeros of a sparse vector. The
	 only purpose of this function is to guard against possible "0"
	 elements delivered by the svec package. This should not happen,
	 so this function should only be seen as a firewall to guard
	 against possible problems if using a sparse-matrix handler that
	 cannot correctly deal with integer matrices.
      */
      template <class T> ivec index_nonzeros(Sparse_Vec<T> x) const;

      //      inline int get_cmax() const {   return (max(sumX1));  } 
      //      inline int get_vmax() const {   return (max(sumX2));  } 
      //      ivec get_coldegree() const; 
      //      ivec get_rowdegree() const; 
    };

  // ---------------------------------------------------------------------------
  // LDPC_Generator_Matrix
  // ---------------------------------------------------------------------------
  /*! \brief LDPC Generator matrix class

  A generator matrix is basically a dense GF(2) matrix with a
  constructor that can create the generator matrix from a parity check
  matrix.

  <b>Please refer to the tutorials for examples of how to use this
  class.</b> 

  \author Erik G. Larsson
  */
  class LDPC_Generator_Matrix 
    { 
      friend class LDPC_Code;
    public:

      //! Default constructor
      LDPC_Generator_Matrix() { type=0; };

      /*!  \brief Construct generator matrix

      This function constructs the generator matrix from a parity check
      matrix (\c LDPC_Parity_Matrix). 

      One can optionally provide an index vector, ind, of variable
      nodes (columns). The algorithm then avoids, to use columns
      corresponding to this index vector. This can be useful for
      example to avoid using variable nodes of a low degree as
      systematic bits.

      \param H parity check matrix
      \param t what type of matrix to generate ("systematic"=systematic 
      form, currently the only supported option)
      \param ind vector of column indices (variable nodes) to avoid in 
      the systematic part

      \note This function modifies the parity check matrix (its
      columns are sorted).
      */
      LDPC_Generator_Matrix(LDPC_Parity_Matrix &H, 
			    std::string t="systematic", ivec ind="");

      //! Assignment operator
      LDPC_Generator_Matrix& operator=(const LDPC_Generator_Matrix &X);
      
    private:
      GF2mat G; // the matrix is stored in transposed form
      int type;  // matrix type (0=not defined, 1=systematic)
      int nvar,ncheck;
    };

  // ---------------------------------------------------------------------------
  // LDPC_Code
  // ---------------------------------------------------------------------------

  /*!
    \brief Low-density parity check (LDPC) codec  

    \author Erik G. Larsson

    This class provides the functionality for encoding and decoding of
    LDPC codes defined via parity check (\c LDPC_Parity_Matrix) and
    generator (\c LDPC_Generator_Matrix) matrices.

    LDPC codecs are constructed from parity and generator
    matrices. Since the procedure of constructing the codec can be
    time-consuming (for example, due to optimization of the parity
    matrix and computation of the generator matrix), codecs can be
    saved to a file for later use.  This class provides functionality
    and a special file format (the file format is designed to optimize
    the operation of the decoder) to do this.

    <b>Please refer to the tutorials \ref ldpc_gen_codes and \ref
    ldpc_bersim_awgn for extensive examples of how to use LDPC
    codes.</b>

    \ingroup fec
  */
  class LDPC_Code : public Channel_Code
    {
    public:
      // ----------------- Constructors ------------------------
      
      //! Default constructor
      LDPC_Code();

      //! Destructor
      virtual ~LDPC_Code() {};

      /*! \brief Constructor, from parity check matrix
    
      \param H the parity check matrix

      This constructor simply calls set_code(). See set_code() for
      documentation.
      */
      LDPC_Code(const LDPC_Parity_Matrix &H);

      /*! \brief Constructor, from parity and generator matrix

      This constructor simply calls set_code(). See set_code() for
      documentation.
      */
      LDPC_Code(const LDPC_Parity_Matrix &H, const LDPC_Generator_Matrix &G);

      /*! \brief Constructor, by reading from file 

      This constructor simply calls set_code(). See set_code() for
      documentation.
      */
      LDPC_Code(std::string filename);

      /*! \brief Set the codec, from parity check matrix
    
      \param H the parity check matrix
      */
      void set_code(const LDPC_Parity_Matrix &H);
  
      /*! \brief Set the codec, from parity and generator matrix

      \param H the parity check matrix
      \param G the generator check matrix
      */
      void set_code(const LDPC_Parity_Matrix &H, const LDPC_Generator_Matrix &G);
      
      /*! \brief Set the codec, by reading from file 

      The file format is defined in the source code and can be read by
      the LDPC_Code::LDPC_Code(string) constructor.  LDPC codecs can
      be saved with the function \c save_to_file.

      \param filename name of the file where the codec is stored
      */
      void set_code(std::string filename);

      /*! \brief Initialize decoder and set parameters

      \param method Currently only supported is "bp" (belief propagation)

      \param options Vector of options to the decoder as follows:

      - options(0): The maximum number of iterations that the decoder 
      performs before terminating. Recommended value 50.

      - options(1): If this flag is 1,  the decoder
      performs a syndrome check each iteration and terminates if valid 
      codeword is found. Recommended value: 1.

      - options(2): If this flag is 1, then the decoder performs a syndrome 
      check before it start decoding. It terminates without
      decoding if LLRin corresponds to a valid codeword (the decoder 
      then sets LLRout=LLRin). Recommended value: 0.

      \param lcalc LLR calculation unit (optional). 

      \note The decoder parameters are not saved to file when saving
      the LDPC code.
      */
      void setup_decoder(std::string method="bp", ivec options="50 1 0",
			 LLR_calc_unit lcalc=LLR_calc_unit()); 

      // ---------------- Save to file ----------------
  
      /*! \brief Save the codec to a file

      \param filename name of the file where to store the codec     
      */
      void save_to_file(std::string filename) const;
  
      // ------------ Encoding  ---------------------
  
      //@{
      /*!  \brief Encode codeword. 

      Currently for systematic generator matrices, this is done by
      brute-force multiplication with the generator matrix.

      \param bitsin vector of input bits
      \param bitsout vector of output bits
      */
      virtual void encode(const bvec &bitsin, bvec &bitsout);
      virtual bvec encode(const bvec &input); 
      //@}

      // ------------ Decoding  ---------------------
      
      //@{
      //! Inherited from base class - not implemented here
      virtual void decode(const bvec &coded_bits, bvec &decoded_bits)
	{ it_error("hard input decoding not implemented"); } ;  
      virtual bvec decode(const bvec &coded_bits)
	{ it_error("hard input decoding not implemented"); return bvec();} ; 
      //@}

      //! This function is a wrapper for bp_decode()
      virtual void decode(const vec &received_signal, bvec &output);

      //! This function is a wrapper for bp_decode()
      virtual bvec decode(const vec &received_signal);

      /*! \brief Belief propagation decoding.  

      This function implements the sum-product message passing decoder
      (Pearl's belief propagation) using LLR values as messages.  A
      fast update mechanism is used for nodes with large degrees.
    
      \param LLRin  vector of input LLR values
      \param LLRout vector of output LLR values

      If the decoder converged to a valid codeword, the function
      returns the number of iterations performed.  Otherwise the
      function returns the number of iterations performed but with
      negative sign.

      Run setup_decoder() first to set number of iterations and other
      parameters. The function uses \c LLR_calc_unit to implement
      table-lookup for the Boxplus operator.  By setting the
      parameters of the \c LLR_calc_unit provided in setup_decoder(),
      one can change the resolution, for example to use a logmax
      approximation. See the documentation of \c LLR_calc_unit for
      details on how to do this.
      */
      int bp_decode(const QLLRvec &LLRin, QLLRvec &LLRout);
   
      /*! \brief Syndrome check, on QLLR vector

      This function performs a syndrome check on a softbit (LLR
      vector). The function returns true for a valid codeword, false
      else.

      \param LLR LLR-vector to check
      */
      bool syndrome_check(const QLLRvec &LLR) const;

      //! Syndrome check, on bit vector
      bool syndrome_check(const bvec &b) const;

      // ------------ Basic information gathering functions ------
  
      //! Get the coderate
      double get_rate(void)  { 
	return   (1.0 - static_cast<double>(ncheck) / nvar); 
      }

      //! Get the number of variable nodes 
      int get_nvar() const { return nvar; }

      //! Get the number of check nodes 
      int get_ncheck() const {return ncheck;}

      //! Get LLR calculation unit used in decoder
      LLR_calc_unit get_llrcalc() const { return llrcalc; }

      //!  Print some properties of the codec in plain text
      friend std::ostream &operator<<(std::ostream &os, LDPC_Code &L);
           
    private:

      // Parity check matrix parameterization 
      int H_is_defined;        // =1 if the parity check matrix defined
      int nvar, ncheck;
      ivec C, V, sumX1, sumX2, iind, jind;  

      // Generator matrix
      LDPC_Generator_Matrix G;  
     
      // decoder parameters
      LLR_calc_unit llrcalc; 
      std::string dec_method;
      ivec dec_options;

      // temporary storage for decoder (memory allocated when codec defined)
      QLLRvec mvc, mcv;

      // Function to compute decoder parameterization 
      void decoder_parameterization(const LDPC_Parity_Matrix &H);
    };

  /*! \relates LDPC_Code
    \brief Print some properties of the LDPC codec in plain text.
  */
  std::ostream &operator<<(std::ostream &os, const LDPC_Code &L);
}

#endif
