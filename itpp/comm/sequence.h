/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2001 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
 * Tobias Ringström, and Jonas Samuelsson.                                   *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*! 
  \file 
  \brief Definitions of binary sequence classes and functions
  \author Tony Ottosson and Pål Frenger

  $Revision$

  $Date$
*/

#ifndef __sequence_h
#define __sequence_h

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/binary.h>

namespace itpp {

  /*! 
    \brief Binary Linear Feedback Shift Register (LFSR)
  
    <ul>
    <li> The LFSR is on Fibonacci form (see e.g. page 104 in Peterson, Ziemer, and Borth
    "Introduction to Spread Spctrum communications", Prentice-Hall, 1995.)</li>
    <li> If the connect_polynomial=1+g1*D+g2*D^2+...+gr*D^r is a primitive polynomial
    a Maximum Length Sequence (m-seq.) of length N=2^r -1 is constructed. Use an
    arbitrary state not equal to zero, to get a phase of the m-sequence.</li>
    <li> For a table of primtive polynomials see e.g. page 117 in the reference above or a suitable book on coding.</li>
    </ul>
  */
  class LFSR {
  public:
    //! Constructor
    LFSR(void) {};
    //! Input connect_polynomial=1+g1*D+g2*D^2+...+gr*D^r in bvec format [g0,g1,...,gr]
    LFSR(const bvec &connections);
    //! Input connect_polynomial=1+g1*D+g2*D^2+...+gr*D^r in octal format
    LFSR(const ivec &connections);
    //! Input connect_polynomial=1+g1*D+g2*D^2+...+gr*D^r in bvec format [g0,g1,...,gr]
    void set_connections(const bvec &connections);
    //! Input connect_polynomial=1+g1*D+g2*D^2+...+gr*D^r in octal format
    void set_connections(const ivec &connections);
    //! Set state (contents in the shift registers) in bvec format
    void set_state(const bvec &state);
    //! Set state (contents in the shift registers) in octal format
    void set_state(const ivec &state);
    //! Shift one step and output binary symbol
    bin shift(void);
    //! Shift no_shifts steps and output bvec
    bvec shift(int no_shifts);
    //! Return length of shift register
    int get_length(void);
    //! Returns the state of the shift register
    bvec get_state(void);
  private:
    bvec memory, Connections;
  };

  /*! 
    \brief Gold Sequences
  */
  class Gold {
  public:
    /*! 
      \brief Class constructor

      Automatic selection of preferred pair of connections. Just give the degree (N=2^degree-1).
      degree=5,7,8,9 available. Only one pair is available for each degree.
    */
    Gold(int degree);
    //! Input connect_polynomials=1+g1*D+g2*D^2+...+gr*D^r in bvec format [g0,g1,...,gr]
    Gold(const bvec &mseq1_connections, const bvec &mseq2_connections);
    //! Input connect_polynomials=1+g1*D+g2*D^2+...+gr*D^r in octal format
    Gold(const ivec &mseq1_connections, const ivec &mseq2_connections);
    //!  Set state (contents in the shift registers) in bvec format
    void set_state(const bvec &state1, const bvec &state2);
    //!  Set state (contents in the shift registers) in octal format
    void set_state(const ivec &state1, const ivec &state2);
    //! Shift one step and output binary symbol
    bin shift(void);
    //! Shift no_shifts steps and output bvec
    bvec shift(int no_shifts);
    //! Returns the length (period) of a Gold-sequence
    int get_sequence_length(void);
    /*! 
      \brief Returns the code family
    
      The Gold code family is defined by the two m-sequences (\a mseq1 and \a mseq2 ) and the sum
      of \a mseq1 and all time shifts of \a mseq2. The return matric thus contain \a N + 2 rows
      and \a N columns, where \a N is the length of the m-sequences.
    */
    bmat get_family(void);
  private:
    int N;
    LFSR mseq1, mseq2;
  };

  // --------------- Inlines ---------------------
  inline bin LFSR::shift(void) {bin temp=memory*Connections;memory.shift_right(temp);return temp;}
  inline int LFSR::get_length(void) {return memory.size();}
  inline bvec LFSR::get_state(void) {return memory;}

  inline bin Gold::shift(void) {return (mseq1.shift()+mseq2.shift());}
  inline int Gold::get_sequence_length(void) {return N;}


  // --------------- Functions ---------------------

  /*! 
    \ingroup miscfunc
    \brief Generates the OVSF (orthogonal variable spreading factor) spreading codes used in WCDMA. 

    The codes are written row-wise in the return matrix.
  */ 
  smat wcdma_spreading_codes(int SF);

} //namespace itpp

#endif // __sequence_h
