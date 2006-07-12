
/*!
 * \file 
 * \brief Class for numerically efficient log-likelihood algebra
 * \author Erik G. Larsson
 *
 * $Date$
 * $Revision $
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


#include <itpp/comm/llr.h>


namespace itpp {

  LLR_calc_unit::LLR_calc_unit()
  {
    init_llr_tables();
  }

  LLR_calc_unit::LLR_calc_unit(short int d1, short int d2, short int d3)
  {
    init_llr_tables(d1,d2,d3);
  }


  void LLR_calc_unit::operator=(const LLR_calc_unit &x)
  {
    Dint1=x.Dint1;
    Dint2=x.Dint2;
    Dint3=x.Dint3;
    logexp_table=x.logexp_table;
  }


  ivec LLR_calc_unit::get_Dint()
  { 
    ivec r(3);
    r(0) = Dint1; 
    r(1) = Dint2;
    r(2) = Dint3;
    return r;
  }

  void LLR_calc_unit::init_llr_tables(short int d1, short int d2, short int d3)
  {
    Dint1 = d1;      // 1<<Dint1 determines how integral LLRs relate to real LLRs (to_double=(1<<Dint)*int_llr)
    Dint2 = d2;      // number of entries in table for LLR operations
    Dint3 = d3;      // table resolution is 2^(-(Dint1-Dint3))
    //    cerr << "Initializing LLR tables, Dint1=" << Dint1 << "   Dint2=" << Dint2 << "  Dint3=" << Dint3  
    //	 << "  resoltion: " << pow(2.0,((double) (Dint3-Dint1))) << endl;
    logexp_table = construct_logexp_table();
  }

  ivec LLR_calc_unit::construct_logexp_table()
  {
    ivec result(Dint2);
    for (int i=0; i<Dint2; i++) {
      double x = pow2(static_cast<double>(Dint3 - Dint1)) * i;
      result(i) = to_qllr(std::log(1 + std::exp(-x)));
    }
    it_assert1(length(result)==Dint2,"Ldpc_codec::construct_logexp_table()");

    return result;
  }

  QLLRvec LLR_calc_unit::to_qllr(const vec &l) {
    int n=length(l);
    ivec result(n);
    for (int i=0; i<n; i++) {
      result.set(i,to_qllr(l(i)));
    }
    return result;
  }
  
  vec LLR_calc_unit::to_double(const QLLRvec &l) {
    int n=length(l);
    vec result(n);
    for (int i=0; i<n; i++) {
      result.set(i,to_double(l(i)));
    }
    return result;
  }

  QLLRmat LLR_calc_unit::to_qllr(const mat &l) {
    int m=l.rows();
    int n=l.cols();
    imat result(m,n);
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
	result.set(i,j,to_qllr(l(i,j)));
      }
    }
    return result;
  }
  
  mat LLR_calc_unit::to_double(const QLLRmat &l) {
    int m=l.rows();
    int n=l.cols();
    mat result(m,n);
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
	result.set(i,j,to_double(l(i,j)));
      }
    }
    return result;
  }
  
  std::ostream &operator<<(std::ostream &os, const LLR_calc_unit &lcu)
  {
    os << "LLR_calc_unit table properties:" << std::endl;
    os << "The granularity in the LLR representation is " 
       << pow2(static_cast<double>(-lcu.Dint1)) << std::endl;
    os << "The LLR scale factor is " << (1 << lcu.Dint1) << std::endl;
    os << "The largest LLR that can be represented is " 
       << lcu.to_double(QLLR_MAX) << std::endl;
    os << "The table resolution is " 
       << pow2(static_cast<double>(lcu.Dint3 - lcu.Dint1)) << std::endl;
    os << "The number of entries in the table is " << lcu.Dint2 << std::endl;
    os << "The tables truncates at the LLR value "
       << pow2(static_cast<double>(lcu.Dint3 - lcu.Dint1)) * lcu.Dint2 
       << std::endl;
    return os;
  }

}
