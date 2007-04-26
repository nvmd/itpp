/*!
 * \file 
 * \brief Implementation of Low-Density Parity Check (LDPC) codes
 * \author Erik G. Larsson, Mattias Andersson and Adam Piatyszek
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
  // LDPC_Parity
  // ---------------------------------------------------------------------------
  
  /*!
    \brief LDPC parity check matrix generic class

    This class provides a basic set of functions needed to represent a
    parity check matrix, which defines an LDPC code. This class is
    used as base class for a set of specific LDPC parity check matrix
    classes, e.g. regular or irregular LDPC codes.

    This class stores a parity check matrix as a sparse matrix.  The
    transpose of the matrix is also stored to enable efficient access
    to its rows.

    All parity check matrices can be loaded from (saved to) a file by
    converting them from (to) a portable \c GF2mat_sparse_alist format.

    However, typically one will want to create a \c LDPC_Code from the
    parity check matrix (and optionally a generator) and save the
    codec binary data instead.

    Please refer to the tutorial \ref ldpc_gen_codes for some examples
    of code generation.
    
    \author Erik G. Larsson, Mattias Andersson and Adam Piatyszek
  */
  class LDPC_Parity {
    friend class LDPC_Code;
  public:     
    //! Default constructor
    LDPC_Parity(): init_flag(false) {}

    //! Constructor that gives an empty matrix of size ncheck x nvar
    LDPC_Parity(int ncheck, int nvar);

    /*! 
      \brief Load an LDPC parity check matrix from a file

      \param filename file name 
      \param format file format

      \note Currently, only "alist" format is supported (see 
      \c GF2mat_sparse_alist for its definition).

      \sa \c load_alist() and \c save_alist()
    */
    LDPC_Parity(const std::string& filename, const std::string& format);

    //! Constructor, from a \c GF2mat_sparse_alist object
    LDPC_Parity(const GF2mat_sparse_alist& alist);

    //! Virtual destructor
    virtual ~LDPC_Parity() {}

    //! Initialize an empty matrix of size ncheck x nvar
    void initialize(int ncheck, int nvar);

    //! Get the parity check matrix, optionally its transposed form
    GF2mat_sparse get_H(bool transpose = false) const {
      return (transpose ? Ht : H);
    } 

    //! Get a specific column from the matrix
    Sparse_Vec<bin> get_col(int c) const { return H.get_col(c); }

    //! Get a specific row from the matrix
    Sparse_Vec<bin> get_row(int r) const { return Ht.get_col(r); }

    //! Get the number of variable nodes (number of columns) 
    int get_nvar() const {
      it_assert_debug(H.cols() == nvar,
		      "LDPC_Parity::get_nvar(): Internal error");
      it_assert_debug(Ht.rows() == nvar,
		      "LDPC_Parity::get_nvar(): Internal error");    
      return nvar; 
    }

    //! Get the number of check nodes (number of rows)
    int get_ncheck() const {
      it_assert_debug(H.rows() == ncheck,
		      "LDPC_Parity::get_ncheck(): Internal error");
      it_assert_debug(Ht.cols() == ncheck,
		      "LDPC_Parity::get_ncheck(): Internal error");
      return ncheck;
    }

    //! Set element (i,j) of the parity check matrix to value
    void set(int i, int j, bin value);

    //! Get element (i,j) of the parity check matrix
    bin get(int i, int j) const {   
      it_assert_debug(H(i,j) == Ht(j,i), "LDPC_Parity::get(): Internal error");
      return H(i,j);  
    }

    //! Get element (i,j) of the parity check matrix
    bin operator()(int i, int j) const {   
      it_assert_debug(H(i,j) == Ht(j,i), 
		      "LDPC_Parity::operator(): Internal error");
      return H(i,j);  
    }

    //! Display some information about the matrix
    virtual void display_stats() const;

    //! Get the code rate
    double get_rate() const { 
      return   (1.0 - static_cast<double>(ncheck) / nvar); 
    }

    //! Import matrix from \c GF2mat_sparse_alist format
    void import_alist(const GF2mat_sparse_alist& H_alist);

    //! Export matrix to \c GF2mat_sparse_alist format
    GF2mat_sparse_alist export_alist() const;

    //! Load matrix from \c alist_file text file in alist format
    void load_alist(const std::string& alist_file);

    //! Save matrix to \c alist_file text file in alist format
    void save_alist(const std::string& alist_file) const;
  
  protected:
    //! Flag that indicates proper initialization
    bool init_flag;
    //! Maximum node degree class can handle
    static const int Nmax = 200;
    //! The parity check matrix
    GF2mat_sparse H;
    //! The transposed parity check matrix
    GF2mat_sparse Ht;
    //! Number of variable nodes
    int nvar; 
    //! Number of check nodes
    int ncheck;
    //! Actual number of ones in each column
    ivec sumX1;
    //! Actual number of ones in each row
    ivec sumX2;
     
    /*! 
      \brief Check for cycles of length L 

      This function implements a recursive routine to find loops. The
      function is mainly a tool for testing and debugging more sophisticated
      functions for graph manipulation.

      \param L length of cycles to look for

      \return The function returns the number of cycles found of length L or
      shorter. Cycles may be counted multiple times.

      \note This function can be very slow for large matrices. It is mainly
      intended as a debugging aid.
    */
    int check_for_cycles(int L) const;

    /*! 
      \brief Check for connectivity between nodes
	
      This function examines whether the point (to_m, to_n) in the matrix
      can be reached from the point (from_m, from_n) using at most L steps.
      A recursive search is used.

      The function can be used to search for cycles in the matrix. To
      search for a cycle of length L, set from_m=to_m and from_n=to_n,
      and godir=0.

      \param from_m starting coordinate, row number
      \param to_m goal coordinate, row number
      \param from_n starting coordinate, column number
      \param to_n goal coordinate, row number
      \param g direction: 1=start going vertically, 2=start 
      going horizontally
      \param L number of permitted steps
     
      \return
      - -1 or -3 : destination unreachable 
      - -2       : meaningless search (started in a "0" point), 
      - -4       : meaningless search
      - >=0      : destination reached with certain number of steps left

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
    int check_connectivity(int from_m, int from_n, int to_m, int to_n,
			   int g, int L) const;

    //      inline int get_cmax() const {   return (max(sumX1));  } 
    //      inline int get_vmax() const {   return (max(sumX2));  } 
    //      ivec get_coldegree() const; 
    //      ivec get_rowdegree() const; 
  };


  // ----------------------------------------------------------------------
  // LDPC_Parity_Unstructured
  // ----------------------------------------------------------------------

  /*!
    \brief Pure abstract class for unstructured LDPC matrices
    
    This class provides a common set of methods for unstructured LDPC
    matrices. For unstructured codes the parity checks are distributed
    at random rather than according to a specific pattern, and the
    generation of a parity matrix can be viewed as drawing a random
    sample from an ensemble of graphs (matrices) that are described by
    a specific degree distribution. 

    This class is used as base for the \c LDPC_Parity_Irregular and \c
    LDPC_Parity_Regular generator classes.

    \author Erik G. Larsson, Mattias Andersson and Adam Piatyszek
  */
  class LDPC_Parity_Unstructured : public LDPC_Parity {
  public:
    //! Display some information about the matrix
    virtual void display_stats() const = 0;

    /*!
      \brief Remove cycles (loops) from unstructured parity check matrix.

      This function implements the cycle removal algorithm presented by
      McGowan and Williamson at the IT workshop 2003. The maximum girth of
      the graph that will be attempted is L. The algorithm is bound to
      remove all loops of length L, insofar this is possible. I.e., it does
      not terminate until it is impossible to remove more cycles by swapping
      two edges.

      \param L Target girth. For example, L=6 attempts to removes all
      4-cycles.
      \return The girth of the graph, i.e. the length of the shortest cycle.
      For example, a return value of 6 means that there are no 4-cycles.

      \note This algorithm can take a long time to run for large L or large
      graphs.
    */
    int cycle_removal_MGW(int L);

  protected:
    //! Generate a random parity check matrix
    void generate_random_H(const ivec& C, const ivec& R, const ivec& cycopt);
  };


  // ----------------------------------------------------------------------
  // LDPC_Parity_Irregular
  // ----------------------------------------------------------------------

  /*!
    \brief Irregular LDPC code generator class
    \author Erik G. Larsson, Mattias Andersson and Adam Piatyszek
  */
  class LDPC_Parity_Irregular : public LDPC_Parity_Unstructured {
  public:
    //! Default constructor
    LDPC_Parity_Irregular() {} 
    //! Constructor that invokes \c generate() method
    LDPC_Parity_Irregular(int Nvar, const vec& var_deg, const vec& chk_deg,
			  const std::string& method = "rand", 
			  const ivec& options = "200 6");

    /*! 
      \brief Generate an irregular LDPC code
       
      \param Nvar number of variable nodes
      \param var_deg vector of variable degree distributions, from an edge
      perspective
      \param chk_deg vector of check degree distributions, from an edge
      perspective
      \param method currently the only provided method is "rand" (see below)
      \param options Determines the level of matrix optimization.  
      
      The "rand" method generates a fully unstructured random
      matrix. Some stochastic optimization is performed to avoid
      cycles. This optimization is controlled via the parameter \c
      options. In particular, the girth for the variable-node-degree-2
      part of the graph can be controlled via the parameter
      options(0). The girth for the rest of the nodes can be
      controlled via the parameter options(1). The recommended value
      for options is "200 6" for graphs of small size, and "100 4" for
      large graphs. The value "0 0" means no optimization. Double
      edges are always avoided.

      \note The "rand" method for graph construction provided in this
      function is not intended to provide the best possible error
      performance, but it is a simple, basic, fast tool that can
      easily be used to build relatively good irregular graphs. Better
      results in terms of performance and error-floor can be achieved
      by using other performance measures for the graph than cycle
      length (for irregular codes, the so-called ACE measure for
      example). Additionally, the "rand" method builds the graph edge
      by edge, with no possibility of removing an edge once it has
      been placed. Better results may be achieved by building the
      graphs by placing pairs or n-tuples of edges at time, for
      example, column by column.

      \note Alternative (user-defined) methods for code generation can
      be implemented by inheriting \c LDPC_Parity_Irregular.
    */
    void generate(int Nvar, const vec& var_deg, const vec& chk_deg, 
		  const std::string& method = "rand", 
		  const ivec& options = "200 6");

    //! Display some information about the matrix
    void display_stats() const { LDPC_Parity::display_stats(); }
  };


  // ----------------------------------------------------------------------
  // LDPC_Parity_Regular
  // ----------------------------------------------------------------------

  /*!
    \brief Regular LDPC code generator class
    \author Erik G. Larsson, Mattias Andersson and Adam Piatyszek
  */
  class LDPC_Parity_Regular : public LDPC_Parity_Unstructured {
  public:
    //! Default constructor
    LDPC_Parity_Regular() {} 
    //! Constructor that invokes \c generate() method
    LDPC_Parity_Regular(int Nvar, int k, int l, 
			const std::string& method = "rand",
			const ivec& options = "200 6");

    /*!
      \brief Generate a (k,l) regular LDPC code 

      \param Nvar  number of variable nodes
      \param k number of ones per column
      \param l number of ones per row
      \param method See \c LDPC_Parity_Irregular::generate()
      \param options See \c LDPC_Parity_Irregular::generate()

      \note Alternative (user-defined) methods for code generation can
      be implemented by inheriting \c LDPC_Parity_Regular.
    */
    void generate(int Nvar, int k, int l, 
		  const std::string& method = "rand",
		  const ivec& options = "200 6");

    //! Display some information about the matrix
    void display_stats() const { LDPC_Parity::display_stats(); }
  };


  // ----------------------------------------------------------------------
  // LDPC_Generator
  // ----------------------------------------------------------------------

  /*!
    \brief LDPC Generator pure virtual base class

    This is an abstract base class for LDPC generators. It provides a
    generic interface that is used by the \c LDPC_Code class. The \c
    LDPC_Generator class can be inherited to create a new type of
    generator. In addition to the default constructor, the following
    three pure virtual methods need to be defined in a derived class:
    \c encode(), \c save() and \c load().

    See the \c LDPC_Generator_Systematic class for an example
    implementation of a derived generator.

    \author Adam Piatyszek
  */
  class LDPC_Generator { 
    friend class LDPC_Code;
  public:
    //! Default constructor
    LDPC_Generator(const std::string& type_in = ""): init_flag(false), 
						     type(type_in) {}
    //! Virtual destructor
    virtual ~LDPC_Generator() {}

    //! Generator specific encode function
    virtual void encode(const bvec &input, bvec &output) = 0;

    //! Return generator type
    std::string get_type() const { return type; }

  protected:
    bool init_flag;		//!< True if generator is initialized
    std::string type;		//!< Generator type

    //! Save generator data to a file
    virtual void save(const std::string& filename) const = 0;
    //! Read generator data from a file
    virtual void load(const std::string& filename) = 0;
  };


  // ----------------------------------------------------------------------
  // LDPC_Generator_Systematic
  // ----------------------------------------------------------------------

  /*!
    \brief Systematic LDPC Generator class

    A generator is basically a dense GF(2) matrix with a constructor
    that can create the generator matrix from a parity check matrix.

    \note Please refer to the tutorials for examples of how to use this
    class.

    \author Erik G. Larsson and Adam Piatyszek
  */
  class LDPC_Generator_Systematic : public LDPC_Generator { 
  public:
    //! Default constructor
    LDPC_Generator_Systematic(): LDPC_Generator("systematic"), G() {}
    //! Parametrized constructor
    LDPC_Generator_Systematic(LDPC_Parity* const H, 
			      bool natural_ordering = false, 
			      const ivec& ind = "");

    //! Virtual destructor
    virtual ~LDPC_Generator_Systematic() {}

    //! Generator specific encode function
    virtual void encode(const bvec &input, bvec &output);

    /*!  
      \brief Construct systematic generator matrix

      This function constructs a systematic generator matrix from a parity
      check matrix (\c LDPC_Parity). The order of the columns is
      randomized unless otherwise requested via the \c natural_ordering
      parameter.

      \param H A pointer to the parity check matrix \c H
      
      \param natural_ordering If this flag is true, the columns are not
      randomly reordered (no interleaving applied), i.e. the function tries
      so far as possible to avoid permuting columns at all. The permutation
      takes place only if absolutely necessary.

      \param ind Vector of column indices (variable nodes) to avoid in the
      systematic part. If this vector is supplied, the algorithm then avoids
      to use variable nodes corresponding to this index vector as systematic
      bits. This can be used for example to avoid using variable nodes of a
      low degree as systematic bits. This parameter is ignored if the
      natural_ordering flag is set.

      \return This function returns the permutation vector \c P on the
      variable nodes that was necessary to construct a full rank generator.
      This is the permutation which effectively has been applied to the
      columns of \c H. The k-th column of the original \c H is the \c P(k)-th
      column of the rearranged \c H.

      \note This function modifies the parity check matrix \c H. Its columns
      may be sorted so that it gets the structure <tt>H = [H1 H2]</tt> where
      \c H2 is square and invertible. The computed generator then satisfies
      <tt>[H1 H2][I; G'] = 0<tt>.
    */
    ivec construct(LDPC_Parity* const H, bool natural_ordering = false, 
		   const ivec& ind = "");

  protected:
    //! Save generator data to a file
    virtual void save(const std::string& filename) const;
    //! Read generator data from a file
    virtual void load(const std::string& filename);
    
  private:
    GF2mat G; // the matrix is stored in transposed form
  };


  // ----------------------------------------------------------------------
  // LDPC_Code
  // ----------------------------------------------------------------------

  /*!
    \ingroup fec
    \brief Low-density parity check (LDPC) codec  

    This class provides the functionality for encoding and decoding of LDPC
    codes defined via \c LDPC_Parity and \c LDPC_Generator classes.

    LDPC codecs are constructed from parity check and generator matrices.
    Since the procedure of constructing the codec can be time-consuming (for
    example, due to optimization of the parity matrix and computation of the
    generator matrix), codecs can be saved to a file for later use. This
    class provides functionality and a special file format (the file format
    is designed to optimize the operation of the decoder) to do this.

    \note Please refer to the tutorials \ref ldpc_gen_codes and \ref
    ldpc_bersim_awgn for extensive examples of how to use LDPC codes.

    \author Erik G. Larsson and Adam Piatyszek
  */
  class LDPC_Code : public Channel_Code {
  public:
    //! Default constructor
    LDPC_Code();

    /*! 
      \brief Constructor, from a parity check matrix and optionally a
      generator
    
      This constructor simply calls \c set_code().
    */
    LDPC_Code(const LDPC_Parity* const H, LDPC_Generator* const G = 0);

    /*!
      \brief Constructor, from a saved file

      This constructor simply calls \c load_code().
    */
    LDPC_Code(const std::string& filename, LDPC_Generator* const G = 0);

    //! Destructor
    virtual ~LDPC_Code() {}

  
    /*!
      \brief Set the codec, from a parity check matrix and optionally a
      generator

      \param H The parity check matrix
      \param G A pointer to the optional generator object
    */
    void set_code(const LDPC_Parity* const H, LDPC_Generator* const G = 0);
      
    /*! 
      \brief Set the codec, by reading from a saved file

      The file format is defined in the source code. LDPC codecs can
      be saved with the function \c save_code().

      \param filename Name of the file where the codec is stored
      \param G A pointer to the optional generator object 

      \note If \c G points at \c 0 (default), the generator date is not
      read from the saved file. This means that the encoding can not be
      performed.
    */
    void load_code(const std::string& filename, LDPC_Generator* const G = 0);

    /*!
      \brief Save the codec to a file

      \param filename Name of the file where to store the codec     

      \note The decoder parameters (\c max_iters, \c syndr_check_each_iter,
      \c syndr_check_at_start and \c llrcalc) are not saved to a file.
    */
    void save_code(const std::string& filename) const;


    /*!
      \brief Set the decoding method
      
      Currently only a belief propagation method ("BP" or "bp") is
      supported.

      \note The default method set in the class constructors is "BP".
    */ 
    void set_decoding_method(const std::string& method);

    /*!
      \brief Set the decoding loop exit conditions
      
      \param max_iters Maximum number of the decoding iterations
      \param stop_if_valid If true, break the decoding loop as soon as
      valid codeword is found. Recommended value: \c true.
      \param syndr_check_at_start If true, perform a syndrome check before
      entering the decoding loop. If LLRin corresponds to a valid codeword,
      set LLRout = LLRin. Recommended value: \c false.

      \note The default values set in the class constructor are: "50",
      "true" and "false", respectively. 
    */
    void set_exit_conditions(int max_iters, 
			     bool syndr_check_each_iter = true,
			     bool syndr_check_at_start = false); 

    //! Set LLR calculation unit
    void set_llrcalc(const LLR_calc_unit& llrcalc);
  

    // ------------ Encoding  ---------------------
  
    /*!
      \brief Encode codeword

      This is a wrapper functions which calls a proper implementation from
      the \c LDPC_Generator object.

      \param input Vector of input bits
      \param output Vector of output bits
    */
    virtual void encode(const bvec &input, bvec &output);
    //! \brief Encode codeword
    virtual bvec encode(const bvec &input); 


    // ------------ Decoding  ---------------------
      
    //! Inherited from the base class - not implemented here
    virtual void decode(const bvec &coded_bits, bvec &decoded_bits)
    { 
      it_error("LDPC_Code::decode(): Hard input decoding not implemented"); 
    }  
    //! Inherited from the base class - not implemented here
    virtual bvec decode(const bvec &coded_bits)
    { 
      it_error("LDPC_Code::decode(): Hard input decoding not implemented"); 
      return bvec();
    } 

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
    double get_rate() const
    { 
      return (1.0 - static_cast<double>(ncheck) / nvar); 
    }

    //! Get the number of variable nodes 
    int get_nvar() const { return nvar; }

    //! Get the number of check nodes 
    int get_ncheck() const { return ncheck; }

    //! Return the decoding method
    std::string get_decoding_method() const { return dec_method; }

    //! Get the maximum number of iterations of the decoder
    int get_nrof_iterations() const { return max_iters; }

    //! Get LLR calculation unit used in decoder
    LLR_calc_unit get_llrcalc() const { return llrcalc; }

    //! Print some properties of the codec in plain text
    friend std::ostream &operator<<(std::ostream &os, const LDPC_Code &C);

  protected:
    bool H_defined;		//!< true if parity check matrix is defined
    bool G_defined;		//!< true if generator is defined
    int nvar;			//!< Number of variable nodes
    int ncheck;			//!< Number of check nodes
    LDPC_Generator *G;		//!< Generator object pointer

    // decoder parameters
    std::string dec_method;	//!< Decoding method 
    int max_iters;		//!< Maximum number of iterations
    bool psc;			//!< check syndrom after each iteration
    bool pisc;			//!< check syndrom before first iteration
    LLR_calc_unit llrcalc;	//!< LLR calculation unit 

    //! Function to compute decoder parameterization 
    void decoder_parameterization(const LDPC_Parity* const H);

    //! Function to check the integrity of the parity check matrix and generator
    void integrity_check();

    //! Initialize decoder
    void setup_decoder(); 

  private:
    // Parity check matrix parameterization 
    ivec C, V, sumX1, sumX2, iind, jind;  

    // temporary storage for decoder (memory allocated when codec defined)
    QLLRvec mvc, mcv;
  };


  /*! 
    \relatesalso LDPC_Code
    \brief Print some properties of the LDPC codec in plain text.
  */
  std::ostream &operator<<(std::ostream &os, const LDPC_Code &C);
}

#endif
