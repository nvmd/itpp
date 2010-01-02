/*!
 * \file
 * \brief Implementation of numerical integration
 * \author Tony Ottosson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include <itpp/base/math/integration.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/help_functions.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/specmat.h>


namespace itpp
{

//! quadstep function
double quadstep(double(*f)(double), double a, double b,
                double fa, double fm, double fb, double is)
{
  double Q, m, h, fml, fmr, i1, i2;
  vec x(2), y(2);

  m = (a + b) / 2;
  h = (b - a) / 4;
  x = vec_2(a + h, b - h);
  y = apply_function<double>(f, x);
  fml = y(0);
  fmr = y(1);

  i1 = h / 1.5 * (fa + 4 * fm + fb);
  i2 = h / 3 * (fa + 4 * (fml + fmr) + 2 * fm + fb);
  i1 = (16 * i2 - i1) / 15;

  if ((is + (i1 - i2) == is) || (m <= a) || (b <= m)) {
    if ((m <= a) || (b <= m)) {
      it_warning("Interval contains no more machine number. Required tolerance may not be met");
    }
    Q = i1;
    return Q;
  }
  else {
    Q = quadstep(f, a, m, fa, fml, fm, is) + quadstep(f, m, b, fm, fmr, fb, is);
  }
  return Q;
}


double quad(double(*f)(double), double a, double b, double tol)
{
  vec x(3), y(3), yy(5);
  double Q, fa, fm, fb, is;

  x = vec_3(a, (a + b) / 2, b);
  y = apply_function<double>(f, x);
  fa = y(0);
  fm = y(1);
  fb = y(2);
  yy = apply_function<double>(f, a + vec(".9501 .2311 .6068 .4860 .8913")
                              * (b - a));
  is = (b - a) / 8 * (sum(y) + sum(yy));

  if (is == 0.0)
    is = b - a;

  is = is * tol / std::numeric_limits<double>::epsilon();
  Q = quadstep(f, a, b, fa, fm, fb, is);

  return Q;
}


//--------------------- quadl() ----------------------------------------

//! quadlstep function
double quadlstep(double(*f)(double), double a, double b,
                 double fa, double fb, double is)
{
  double Q, h, m, alpha, beta, mll, ml, mr, mrr, fmll, fml, fm, fmr, fmrr,
  i1, i2;
  vec x(5), y(5);

  h = (b - a) / 2;
  m = (a + b) / 2;
  alpha = std::sqrt(2.0 / 3);
  beta = 1.0 / std::sqrt(5.0);
  mll = m - alpha * h;
  ml = m - beta * h;
  mr = m + beta * h;
  mrr = m + alpha * h;
  x(0) = mll;
  x(1) = ml;
  x(2) = m;
  x(3) = mr;
  x(4) = mrr;

  y = apply_function<double>(f, x);

  fmll = y(0);
  fml = y(1);
  fm = y(2);
  fmr = y(3);
  fmrr = y(4);

  i2 = (h / 6) * (fa + fb + 5 * (fml + fmr));
  i1 = (h / 1470) * (77 * (fa + fb) + 432 * (fmll + fmrr) + 625 * (fml + fmr) + 672 * fm);

  if ((is + (i1 - i2) == is) || (mll <= a) || (b <= mrr)) {
    if ((m <= a) || (b <= m)) {
      it_warning("Interval contains no more machine number. Required tolerance may not be met");
    }
    Q = i1;
    return Q;
  }
  else {
    Q = quadlstep(f, a, mll, fa, fmll, is) + quadlstep(f, mll, ml, fmll, fml, is) + quadlstep(f, ml, m, fml, fm, is) +
        quadlstep(f, m, mr, fm, fmr, is) + quadlstep(f, mr, mrr, fmr, fmrr, is) + quadlstep(f, mrr, b, fmrr, fb, is);
  }
  return Q;
}

double quadl(double(*f)(double), double a, double b, double tol)
{
  double Q, m, h, alpha, beta, x1, x2, x3, fa, fb, i1, i2, is, s, erri1, erri2, R;
  vec x(13), y(13);
  double tol2 = tol;

  m = (a + b) / 2;
  h = (b - a) / 2;

  alpha = std::sqrt(2.0 / 3);
  beta = 1.0 / std::sqrt(5.0);

  x1 = .942882415695480;
  x2 = .641853342345781;
  x3 = .236383199662150;
  x(0) = a;
  x(1) = m - x1 * h;
  x(2) = m - alpha * h;
  x(3) = m - x2 * h;
  x(4) = m - beta * h;
  x(5) = m - x3 * h;
  x(6) = m;
  x(7) = m + x3 * h;
  x(8) = m + beta * h;
  x(9) = m + x2 * h;
  x(10) = m + alpha * h;
  x(11) = m + x1 * h;
  x(12) = b;

  y = apply_function<double>(f, x);

  fa = y(0);
  fb = y(12);
  i2 = (h / 6) * (y(0) + y(12) + 5 * (y(4) + y(8)));
  i1 = (h / 1470) * (77 * (y(0) + y(12)) + 432 * (y(2) + y(10)) + 625 * (y(4) + y(8)) + 672 * y(6));

  is = h * (.0158271919734802 * (y(0) + y(12)) + .0942738402188500 * (y(1) + y(11)) + .155071987336585 * (y(2) + y(10)) +
            .188821573960182 * (y(3) + y(9)) + .199773405226859 * (y(4) + y(8)) + .224926465333340 * (y(5) + y(7)) + .242611071901408 * y(6));

  s = sign(is);
  if (s == 0.0)
    s = 1;

  erri1 = std::abs(i1 - is);
  erri2 = std::abs(i2 - is);

  R = 1;
  if (erri2 != 0.0)
    R = erri1 / erri2;

  if (R > 0 && R < 1)
    tol2 = tol2 / R;

  is = s * std::abs(is) * tol2 / std::numeric_limits<double>::epsilon();
  if (is == 0.0)
    is = b - a;

  Q = quadlstep(f, a, b, fa, fb, is);

  return Q;
}


} // namespace itpp
