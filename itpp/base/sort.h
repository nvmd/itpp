/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2004 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Sorting functions
  \author Tony Ottosson

  1.21

  2003/06/17 11:48:36
*/

#ifndef __sort_h
#define __sort_h

#include "itpp/base/vec.h"
//#include "base/mat.h"
//#include "base/converters.h"
//#include "base/scalfunc.h"
#include "itpp/base/itassert.h"
//#include "base/specmat.h"
//#include "base/binary.h"

namespace itpp {


  /*! 
    \relates Vec
    \brief Sort the the vector in increasing order
  */
  template<class T>
    void sort(Vec<T> &data);

  /*! 
    \relates Vec
    \brief Return an index vector corresponding to a sorted vector (increasing order)
  */
  template<class T>
    ivec sort_index(const Vec<T> &data);



  template<class T> void QS(int low, int high, Vec<T> &data) {
    int plow, phigh;
    T a,test;

    if (high>low) {
      a=data[low];
      plow=low;
      phigh=high;
      test=data[phigh];
      while (plow<phigh) {
	if (test<a) {
	  data[plow]=test;
	  plow++;
	  test=data[plow];
	} else {
	  data[phigh]=test;
	  phigh--;
	  test=data[phigh];
	}
      }
      data[plow]=a;
      QS(low,plow-1,data);
      QS(plow+1,high,data);
    }
  }

  template<class T>
    void sort(Vec<T> &data)
    {
      QS(0,data.size()-1,data);
    }

  template<class T>
    void QSindex(int low, int high, ivec &indexlist, const Vec<T> &data)
    {
      int plow,phigh,testindex,aindex;
      T a,test;

      if (high>low) {
	aindex=indexlist[low];
	a=data[aindex];
	plow=low;
	phigh=high;
	testindex=indexlist[phigh];
	test=data[testindex];
	while (plow<phigh) {
	  if (test<a) {
	    indexlist[plow]=testindex;
	    plow++;
	    testindex=indexlist[plow];
	    test=data[testindex];
	  } else {
	    indexlist[phigh]=testindex;
	    phigh--;
	    testindex=indexlist[phigh];
	    test=data[testindex];
	  }
	}
	indexlist[plow]=aindex;
	QSindex(low,plow-1,indexlist,data);
	QSindex(plow+1,high,indexlist,data);
      }
    }

  template<class T>
    ivec sort_index(const Vec<T> &data)
    {
      int N=data.length(),i;
      ivec indexlist(N);

      for(i=0;i<N;i++) {
	indexlist(i)=i;
      }
      QSindex(0,N-1,indexlist,data);
      return indexlist;
    }


} //namespace itpp

#endif // __sort_h

