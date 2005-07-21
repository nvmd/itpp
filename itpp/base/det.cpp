/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2005 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Implementation of determinant calculations
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#include <itpp/base/vec.h>
#include <itpp/base/itassert.h>
#include <itpp/base/lu.h>

namespace itpp {


  /* Determinant of square matrix.
     Calculate determinant of inmatrix (Uses LU-factorisation)
     (See Theorem 3.2.1 p. 97 in Golub & van Loan, "Matrix Computations").

     det(X) = det(P')*det(L)*det(U) = det(P')*1*prod(diag(U))
  */
  double det(const mat &X)
  {
    it_assert1(X.rows()==X.cols(),"det : Only square matrices");

    mat L, U;
    ivec p;
    double s=1.0;
    int i;

    lu(X,L,U,p); // calculate LU-factorisation

    double temp=U(0,0);
    for (i=1;i<X.rows();i++) {
      temp*=U(i,i);
    }

    // Calculate det(P´). Equal to (-1)^(no row changes)
    for (i=0; i<p.size(); i++)
      if (i != p(i))
	s *=-1.0;

    return temp*s;
  }


  /* Determinant of complex square matrix.
     Calculate determinant of inmatrix (Uses LU-factorisation)
     (See Theorem 3.2.1 p. 97 in Golub & van Loan, "Matrix Computations").

     det(X) = det(P')*det(L)*det(U) = det(P')*1*prod(diag(U))

     Needs LU-factorization of complex matrices (LAPACK)
  */
  std::complex<double> det(const cmat &X)
  {
    it_assert1(X.rows()==X.cols(),"det : Only square matrices");

    int i;
    cmat L, U;
    ivec p;
    double s=1.0;

    lu(X,L,U,p); // calculate LU-factorisation

    std::complex<double> temp=U(0,0);
    for (i=1;i<X.rows();i++) {
      temp*=U(i,i);
    }

    // Calculate det(P´). Equal to (-1)^(no row changes)
    for (i=0; i<p.size(); i++)
      if (i != p(i))
	s *=-1.0;

    return temp*s;
  }


} //namespace itpp
