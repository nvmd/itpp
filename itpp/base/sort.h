/*!
 * \file
 * \brief Sorting functions
 * \author Tony Ottosson
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

#ifndef SORT_H
#define SORT_H

#include <itpp/base/vec.h>


namespace itpp {

  /*! 
    \relates Vec
    \brief Sort the \c data vector in increasing order
  */
  template<class T>
  void sort(Vec<T> &data)
  {
    QS(0,data.size()-1,data);
  }

  /*! 
    \relates Vec
    \brief Return an index vector corresponding to a sorted vector 
    \c data (increasing order)
  */
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

  /*! 
    \relates Vec
    \brief Quick sort function of a subset of a vector \c data 

    \param low Start index of a subvector to be sorted
    \param high End index of a subvector to be sorted
    \param data Data vector, in which a part of it is to be sorted

    \note This function uses recurence.
  */
  template<class T>
  void QS(int low, int high, Vec<T> &data)
  {
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

  /*! 
    \relates Vec
    \brief Quick sort function, which gives a sorted index vector
    \c indexlist 

    \param low Start index of a subvector to be sorted
    \param high End index of a subvector to be sorted
    \param data Data vector, in which a part of it is to be sorted
    \param indexlist Result of sorting in the form of sorted indexes

    \note This function uses recurence.
  */
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

} // namespace itpp

#endif // #ifndef SORT_H

