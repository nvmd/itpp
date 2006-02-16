
/*!
 * \file 
 * \brief Class for numerically efficient log-likelihood algebra
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


#ifndef LLR_H
#define LLR_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/specmat.h>
#include <itpp/base/matfunc.h>


namespace itpp {
  
  /*! \relates LLR_calc_unit 
    The quantized LLR (QLLR) representation, scalar form 
  */
  typedef signed int QLLR;

  /*!  \relates LLR_calc_unit 
    The quantized LLR (QLLR) representation, vector form
  */
  typedef Vec<QLLR> QLLRvec;

  /*!  \relates LLR_calc_unit 
    The quantized LLR (QLLR) representation, matrix form
  */
  typedef Mat<QLLR> QLLRmat;
  
  /*!  \relates LLR_calc_unit 
    The largest possible QLLR value
  */
  const QLLR QLLR_MAX=(INT_MAX>>2) ;
 // added some margin to make sure the sum of two LLR is still permissible

  /*! 
    \brief Log-likelihood algebra calculation unit.

    This class contains functions for algebra with log-likelihood
    ratios (LLRs). The (sole) purpose of this class is to provide
    numerically efficient functions for turbo and LDPC codes, which
    rely on certain nonlinear operations on LLRs.

    An LLR for an information bit b is defined according to \f[ L =
    \log \frac{P(b=1)}{P(b=0)} \f] and it is in general a real number.  In the
    class, LLR values are represented via the special type, "quantized
    LLR" (QLLR).  The relation between the quantized representation
    and the real (floating-point) LLR value is \f[ \mbox{QLLR} = \mbox{round}
    ((1<<\mbox{Dint1})*\mbox{LLR}) \f]  The user parameter Dint1 determines the
    granularity of the quantization, and it can be set arbitrarily.
    The functions to_double() and to_qllr() can be used to perform
    conversions between the two representations (QLLR to
    floating-point, and vice versa).

    The class provides functions for the computation of the Jacobian
    logarithm and Hagenauer's "boxplus" operator.  These functions are
    based on a table-lookup.  The resolution of the table is
    determined by the parameters Dint2 and Dint3. See the class
    constructor for more detail.  When an object of LLR_calc_unit is
    created, corresponding lookup-tables are also generated. The
    resolution of these tables can be adjusted by providing parameters
    to the constructor.

    The variable table resolution allows one to study complexity
    versus accuracy (i.e., how different table resolutions would
    degrade performance) to some extent. Yet the main purpose of the
    QLLR representation is to provide a tool for writing efficient
    simulation code, rather than to provide for bit-level
    (fixed-point) simulations. For bit-level simulations, a true fixed
    point representation of LLRs would be preferred/required.  With
    the default setting of the table parameters, using the QLLR type
    is practically as accurate (but much faster) as using "double" to
    represent LLRs.  Decoder implementations may then provide
    functions using QLLR, fixed-point, or double (for compatibility
    reasons) representations of LLR values.

    Note: the QLLR type does not check that the correct quantization
    level is used. I.e., in theory it would be possible to add two
    QLLR types with different quantization (Dint) parameters.  This is
    intentionally implemented this way to achieve maximum runtime
    efficiency.

  */
  class LLR_calc_unit {
  public:
    //! Constructor, using the default table resolution
    LLR_calc_unit();

    /*! \brief Constructor, using a specific table resolution. 

    See init_llr_tables() for more detail on the parameters.
     */
    LLR_calc_unit(const short int in_Dint1, const short int in_Dint2, const short int in_Dint3);

    /*! \brief Set the quantization and table parameters

      \param Dint1 Determines the relation between LLR represented as
      real number and as integer.  The relation is \f[ \mbox{QLLR} = \mbox{round}
      ((1<<\mbox{Dint1})*\mbox{LLR}) \f]

      \param Dint2 Number of entries in the table. If this is zero,
      then logmap becomes logmax.

      \param Dint3 Determines the table resolution. The spacing between each 
      entry is \f[ 2^{-(Dint1-Dint3)} \f]

      The default parameter values are chosen to give a performance
      practically indistinguishable from that of using floating point
      calculations.
    */
    void init_llr_tables(const short int d1=12, const short int d2=300, const short int d3=7);
    //    void init_llr_tables(const short int d1=10, const short int d2=300, const short int d3=5);

    //! Convert a "real" LLR value to an LLR type
    inline QLLR to_qllr(const double &l);

    //! Convert a vector of "real" LLR values to an LLR type
    QLLRvec to_qllr(const vec &l);

    //! Convert a matrix of "real" LLR values to an LLR type
    QLLRmat to_qllr(const mat &l);
  
    //! Convert an LLR type to a "real" LLR 
    inline double to_double(const QLLR &l) const { return ( ((double) l) / ((double) (1<<Dint1))); };

    //! Convert a vector of LLR types to a "real" LLR 
    vec to_double(const QLLRvec &l);

    //! Convert a matrix of LLR types to a "real" LLR 
    mat to_double(const QLLRmat &l);

    /*! \brief Jacobian logarithm. 

    This function computes \f[ \log(\exp(a)+\exp(b)) \f]
    */
    inline QLLR jaclog(QLLR a, QLLR b);
    // Note: a version of this function taking "double" values as input
    // is deliberately omitted, because this is rather slow. 

    /*! \brief Hagenauer's "Boxplus" operator.  

      This function computes
      \f[ \mbox{sign}(a)*\mbox{sign}(b)*\mbox{min}(|a|,|b|)+f(|a+b|)-f(|a-b|) \f]
      where \f[ f(x) = \log(1+\exp(-x))  \f]
    */
    inline QLLR Boxplus(QLLR a, QLLR b);

    /*! \brief Logexp operator.

    This function computes \f[ f(x) = \log(1+\exp(-x)) \f]  */
    QLLR logexp(QLLR x);

    //! Retrieve the table resolution values
    ivec get_Dint();
  
    //! Assignment operator for LLR_calc_unit 
    void operator=(const LLR_calc_unit&);

    //!  Print some properties of the LLR calculation unit (i.e., the lookup table parameters) in plain text
    friend std::ostream &operator<<(std::ostream &os, const LLR_calc_unit &lcu);
    
  private:
    //! Compute the table for \f[ f(x) = \log(1+\exp(-x)) \f]
    ivec construct_logexp_table();
    
    //! The lookup tables for the decoder
    ivec logexp_table;
    
    //! Decoder (lookup-table) parameters
    short int Dint1, Dint2, Dint3;

  };

  /*! \relates LLR_calc_unit
    \brief Print some properties of the LLR calculation unit in plain text.
  */
  std::ostream &operator<<(std::ostream &os, const LLR_calc_unit &lcu);

  

  // ================ IMPLEMENTATIONS OF SOME LIKELIHOOD ALGEBRA FUNCTIONS =========== 

  inline QLLR LLR_calc_unit::to_qllr(const double &l)
    { 
      it_assert0(l<=to_double(QLLR_MAX),"LLR_calc_unit::to_qllr(): overflow");
      it_assert0(l>=-to_double(QLLR_MAX),"LLR_calc_unit::to_qllr(): overflow");
      return ( (int) std::floor(0.5+(1<<Dint1)*l) ); 
    };
  
  inline QLLR LLR_calc_unit::Boxplus(QLLR a, QLLR b)
    {
      /* These boundary checks are meant to completely eliminate the
	 possibility of any numerical problem.  Perhaps they are not
	 strictly necessary.  */
      if (a>=QLLR_MAX) { 
	return b;
      }
      if (b>=QLLR_MAX) {
	return a;
      }
      if (a<=-QLLR_MAX) {
	return (-b);
      }
      if (b<=-QLLR_MAX) {
	return (-a);
      }

      QLLR a_abs = (a>0 ? a : -a);
      QLLR b_abs = (b>0 ? b : -b);
      QLLR minabs = (a_abs>b_abs ? b_abs : a_abs);
      QLLR term1 = (a>0 ? (b>0 ? minabs : -minabs) : (b>0 ? -minabs : minabs));
      QLLR apb = a+b;
      QLLR term2 = logexp((apb>0 ? apb : -apb));
      QLLR amb = a-b;
      QLLR term3 = logexp((amb>0 ? amb : -amb));

      QLLR result = term1+term2-term3;
      if (result>=QLLR_MAX) {         
	return QLLR_MAX; 
      } else if (result<=-QLLR_MAX) {
	return (-QLLR_MAX);
      } else {
	return result;
      }
    }

  inline QLLR LLR_calc_unit::logexp(QLLR x)
    {
      it_assert0(x>=0,"LLR_calc_unit::logexp() is not defined for negative LLR values");
      int ind = x>>Dint3;
      if (ind>=Dint2) {
	return 0;
      }

      it_assert0(ind>=0,"LLR_calc_unit::logexp() internal error");
      it_assert0(ind<Dint2,"LLR_calc_unit::logexp() internal error");

      // With interpolation 
      // int delta=x-(ind<<Dint3);
      // return ((delta*logexp_table(ind+1) + ((1<<Dint3)-delta)*logexp_table(ind)) >> Dint3);

      // Without interpolation
      return logexp_table(ind);
    }

  inline QLLR LLR_calc_unit::jaclog(QLLR a, QLLR b )
    {
      QLLR x,maxab;
  
      if (a>b) {
	maxab = a;
	x=a-b;
      } else {
	maxab = b;
	x=b-a;
      }
 
      if (maxab>=QLLR_MAX) {
	return QLLR_MAX; 
      } else {
	return (maxab + logexp(x));
      }
    };

}

#endif
