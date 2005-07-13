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
  \brief Implementation of special operations on vectors and matricies optimized for speed
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#include "itpp/base/binary.h"
#include "itpp/base/fastmath.h"
 
namespace itpp { 

  // m=m-v*v'*m
  void sub_v_vT_m(mat &m, const vec &v)
  {
    vec v2(m.cols());
    double tmp, *v2p;
    const double *vp;
    int i, j;

    it_assert(v.size() == m.rows(), "sub_v_vT_m()");

    v2p = v2._data();
    for (j=0; j<m.cols(); j++) {
      tmp = 0.0;
      vp=v._data();
      for (i=0; i<m.rows(); i++)
	tmp += *(vp++) * m._elem(i,j);
      *(v2p++) = tmp;
    }

    vp=v._data();
    for (i=0; i<m.rows(); i++) {
      v2p = v2._data();
      for (j=0; j<m.cols(); j++)
	m._elem(i,j) -= *vp * *(v2p++);
      vp++;
    }
  }

  // m=m-m*v*v'
  void sub_m_v_vT(mat &m, const vec &v)
  {
    vec v2(m.rows());
    double tmp, *v2p;
    const double *vp;
    int i, j;

    it_assert(v.size() == m.cols(), "sub_m_v_vT()");
    
    v2p = v2._data();
    for (i=0; i<m.rows(); i++) {
      tmp = 0.0;
      vp = v._data();
      for (j=0; j<m.cols(); j++)
	tmp += *(vp++) * m._elem(i,j);
      *(v2p++) = tmp;
    }

    v2p = v2._data();
    for (i=0; i<m.rows(); i++) {
      vp=v._data();
      for (j=0; j<m.cols(); j++)
	m._elem(i,j) -= *v2p * *(vp++);
      v2p++;
    }
  }

} //namespace itpp
