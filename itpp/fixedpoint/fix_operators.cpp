/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 2005 by Johan Bergman.                                      *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*!
  \file
  \brief Operators for Fix, Fixed, CFix and CFixed
  \author Johan Bergman
  
  $Revision$
  
  $Date$
*/

#include <itpp/fixedpoint/fix_operators.h>
#include <itpp/base/itassert.h>

namespace itpp {

  /////////////////////////////////
  // Operators for Fix and Fixed //
  /////////////////////////////////

  Fix operator+(const Fix &x, const Fix &y)
  {
    return Fix(x.get_re() + y.get_re(),
               assert_shifts(x, y),
               0, 0);
  }

  Fix operator-(const Fix &x, const Fix &y)
  {
    return Fix(x.get_re() - y.get_re(),
               assert_shifts(x, y),
               0, 0);
  }

  Fix operator*(const Fix &x, const Fix &y)
  {
    return Fix(x.get_re() * y.get_re(),
               x.get_shift() + y.get_shift(),
               0, 0);
  }

  Fix operator/(const Fix &x, const Fix &y)
  {
    return Fix(x.get_re() / y.get_re(),
               x.get_shift() - y.get_shift(),
               0, 0);
  }

  Fix operator+(const Fix &x, const int y)
  {
    return Fix(x.get_re() + y,
               assert_shifts(x, y),
               0, 0);
  }

  Fix operator-(const Fix &x, const int y)
  {
    return Fix(x.get_re() - y,
               assert_shifts(x, y),
               0, 0);
  }

  Fix operator*(const Fix &x, const int y)
  {
    return Fix(x.get_re() * y,
               x.get_shift(),
               0, 0);
  }

  Fix operator/(const Fix &x, const int y)
  {
    return Fix(x.get_re() / y,
               x.get_shift(),
               0, 0);
  }

  Fix operator+(const int x, const Fix &y)
  {
    return Fix(x + y.get_re(),
               assert_shifts(y, x),
               0, 0);
  }

  Fix operator-(const int x, const Fix &y)
  {
    return Fix(x - y.get_re(),
               assert_shifts(y, x),
               0, 0);
  }

  Fix operator*(const int x, const Fix &y)
  {
    return Fix(x * y.get_re(),
               y.get_shift(),
               0, 0);
  }

  Fix operator/(const int x, const Fix &y)
  {
    return Fix(x / y.get_re(),
               -y.get_shift(),
               0, 0);
  }

  //#ifdef _MSC_VER
  // These operators are not defined when using the Microsoft Developer Studio compiler due
  // to internal compiling error. Hopefully this compiler-bug will be fixed in the near future.
  //#else

  fixvec operator+(const fixvec &a, const ivec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes do not match");
    fixvec temp(a);
    for (int i=0; i<a.size(); i++) {
      temp(i) += b(i);
    }
    return temp;
  }

  Fix operator*(const fixvec &a, const ivec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes do not match");
    Fix temp(0);
    for (int i=0; i<a.size(); i++) {
      temp += a(i) * b(i);
    }
    return temp;
  }

  fixmat operator+(const fixmat &a, const imat &b)
  {
    it_assert1(a.cols()==b.cols() && a.rows()==b.rows(), "operator+(): sizes do not match");
    fixmat temp(a);

    for (int i=0; i<a.rows(); i++) {
      for (int j=0; j<a.cols(); j++) {
        temp(i,j) += b(i,j);
      }
    }
    return temp;
  }

  fixmat operator*(const fixmat &a, const imat &b)
  {
    it_assert1(a.cols() == b.rows(), "operator*: wrong sizes");
    fixmat r(a.rows(), b.cols());
  
    Fix tmp;
    int i, j, k;
    Fix *tr=r._data();
    const Fix *t1;
    const int *t2=b._data();
  
    for (i=0; i<r.cols(); i++) {
      for (j=0; j<r.rows(); j++) {
        tmp = Fix(0); t1 = a._data()+j;
        for (k=a.cols(); k>0; k--) {
          tmp += *(t1) * *(t2++);
          t1 += a.rows();
        }
        *(tr++) = tmp; t2 -= b.rows();
      }
      t2 += b.rows();
    }  
    return r;
  }

  //#endif // _MSC_VER

  ///////////////////////////////////
  // Operators for CFix and CFixed //
  ///////////////////////////////////

  CFix operator+(const CFix &x, const CFix &y)
  {
    return CFix(x.get_re() + y.get_re(),
                x.get_im() + y.get_im(),
                assert_shifts(x, y),
                0, 0);
  }

  CFix operator-(const CFix &x, const CFix &y)
  {
    return CFix(x.get_re() - y.get_re(),
                x.get_im() - y.get_im(),
                assert_shifts(x, y),
                0, 0);
  }

  CFix operator*(const CFix &x, const CFix &y)
  {
    return CFix(x.get_re()*y.get_re() - x.get_im()*y.get_im(),
                x.get_re()*y.get_im() + x.get_im()*y.get_re(),
                x.get_shift() + y.get_shift(),
                0, 0);
  }

  CFix operator/(const CFix &x, const CFix &y)
  {
    fixrep denominator = y.get_re()*y.get_re() + y.get_im()*y.get_im();
    return CFix((x.get_re()*y.get_re() + x.get_im()*y.get_im())/denominator,
                (x.get_im()*y.get_re() - x.get_re()*y.get_im())/denominator,
                x.get_shift() - y.get_shift(),
                0, 0);
  }

  CFix operator+(const CFix &x, const Fix &y)
  {
    return CFix(x.get_re() + y.get_re(),
                x.get_im(),
                assert_shifts(x, y),
                0, 0);
  }

  CFix operator-(const CFix &x, const Fix &y)
  {
    return CFix(x.get_re() - y.get_re(),
                x.get_im(),
                assert_shifts(x, y),
                0, 0);
  }

  CFix operator*(const CFix &x, const Fix &y)
  {
    return CFix(x.get_re() * y.get_re(),
                x.get_im() * y.get_re(),
                x.get_shift() + y.get_shift(),
                0, 0);
  }

  CFix operator/(const CFix &x, const Fix &y)
  {
    return CFix(x.get_re() / y.get_re(),
                x.get_im() / y.get_re(),
                x.get_shift() - y.get_shift(),
                0, 0);
  }

  CFix operator+(const Fix &x, const CFix &y)
  {
    return CFix(x.get_re() + y.get_re(),
                y.get_im(),
                assert_shifts(y, x),
                0, 0);
  }

  CFix operator-(const Fix &x, const CFix &y)
  {
    return CFix(x.get_re() - y.get_re(),
                -y.get_im(),
                assert_shifts(y, x),
                0, 0);
  }

  CFix operator*(const Fix &x, const CFix &y)
  {
    return CFix(x.get_re() * y.get_re(),
                x.get_re() * y.get_im(),
                x.get_shift() + y.get_shift(),
                0, 0);
  }

  CFix operator/(const Fix &x, const CFix &y)
  {
    fixrep denominator = y.get_re()*y.get_re() + y.get_im()*y.get_im();
    return CFix(x.get_re() * y.get_re() / denominator,
                -x.get_re() * y.get_im() / denominator,
                x.get_shift() - y.get_shift(),
                0, 0);
  }

  CFix operator+(const CFix &x, const int y)
  {
    return CFix(x.get_re() + y,
                x.get_im(),
                assert_shifts(x, y),
                0, 0);
  }

  CFix operator-(const CFix &x, const int y)
  {
    return CFix(x.get_re() - y,
                x.get_im(),
                assert_shifts(x, y),
                0, 0);
  }

  CFix operator*(const CFix &x, const int y)
  {
    return CFix(x.get_re() * y,
                x.get_im() * y,
                x.get_shift(),
                0, 0);
  }

  CFix operator/(const CFix &x, const int y)
  {
    return CFix(x.get_re() / y,
                x.get_im() / y,
                x.get_shift(),
                0, 0);
  }

  CFix operator+(const int x, const CFix &y)
  {
    return CFix(x + y.get_re(),
                y.get_im(),
                assert_shifts(y, x),
                0, 0);
  }

  CFix operator-(const int x, const CFix &y)
  {
    return CFix(x - y.get_re(),
                -y.get_im(),
                assert_shifts(y, x),
                0, 0);
  }

  CFix operator*(const int x, const CFix &y)
  {
    return CFix(x * y.get_re(),
                x * y.get_im(),
                y.get_shift(),
                0, 0);
  }

  CFix operator/(const int x, const CFix &y)
  {
    fixrep denominator = y.get_re()*y.get_re() + y.get_im()*y.get_im();
    return CFix(x * y.get_re() / denominator,
                -x * y.get_im() / denominator,
                -y.get_shift(),
                0, 0);
  }

  //#ifdef _MSC_VER
  // These operators are not defined when using the Microsoft Developer Studio compiler due
  // to internal compiling error. Hopefully this compiler-bug will be fixed in the near future.
  //#else

  cfixvec operator+(const cfixvec &a, const fixvec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes do not match");
    cfixvec temp(a);
    for (int i=0; i<a.size(); i++) {
      temp(i) += b(i);
    }
    return temp;
  }

  CFix operator*(const cfixvec &a, const fixvec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes do not match");
    CFix temp(0);
    for (int i=0; i<a.size(); i++) {
      temp += a(i) * b(i);
    }
    return temp;
  }

  cfixmat operator+(const cfixmat &a, const fixmat &b)
  {
    it_assert1(a.cols()==b.cols() && a.rows()==b.rows(), "operator+(): sizes do not match");
    cfixmat temp(a);

    for (int i=0; i<a.rows(); i++) {
      for (int j=0; j<a.cols(); j++) {
        temp(i,j) += b(i,j);
      }
    }
    return temp;
  }

  cfixmat operator*(const cfixmat &a, const fixmat &b)
  {
    it_assert1(a.cols() == b.rows(), "operator*: wrong sizes");
    cfixmat r(a.rows(), b.cols());
  
    CFix tmp;
    int i, j, k;
    CFix *tr=r._data();
    const CFix *t1;
    const Fix *t2=b._data();
  
    for (i=0; i<r.cols(); i++) {
      for (j=0; j<r.rows(); j++) {
        tmp = CFix(0); t1 = a._data()+j;
        for (k=a.cols(); k>0; k--) {
          tmp += *(t1) * *(t2++);
          t1 += a.rows();
        }
        *(tr++) = tmp; t2 -= b.rows();
      }
      t2 += b.rows();
    }  
    return r;
  }

  cfixvec operator+(const cfixvec &a, const ivec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes do not match");
    cfixvec temp(a);
    for (int i=0; i<a.size(); i++) {
      temp(i) += b(i);
    }
    return temp;
  }

  CFix operator*(const cfixvec &a, const ivec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes do not match");
    CFix temp(0);
    for (int i=0; i<a.size(); i++) {
      temp += a(i) * b(i);
    }
    return temp;
  }

  cfixmat operator+(const cfixmat &a, const imat &b)
  {
    it_assert1(a.cols()==b.cols() && a.rows()==b.rows(), "operator+(): sizes do not match");
    cfixmat temp(a);

    for (int i=0; i<a.rows(); i++) {
      for (int j=0; j<a.cols(); j++) {
        temp(i,j) += b(i,j);
      }
    }
    return temp;
  }

  cfixmat operator*(const cfixmat &a, const imat &b)
  {
    it_assert1(a.cols() == b.rows(), "operator*: wrong sizes");
    cfixmat r(a.rows(), b.cols());
  
    CFix tmp;
    int i, j, k;
    CFix *tr=r._data();
    const CFix *t1;
    const int *t2=b._data();
  
    for (i=0; i<r.cols(); i++) {
      for (j=0; j<r.rows(); j++) {
        tmp = CFix(0); t1 = a._data()+j;
        for (k=a.cols(); k>0; k--) {
          tmp += *(t1) * *(t2++);
          t1 += a.rows();
        }
        *(tr++) = tmp; t2 -= b.rows();
      }
      t2 += b.rows();
    }  
    return r;
  }

  //#endif // _MSC_VER

} //namespace itpp
