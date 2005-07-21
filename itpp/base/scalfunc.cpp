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
  \brief Implementation of scalar functions
  \author Tony Ottosson and Pål Frenger

  1.18

  2003/05/22 08:55:19
*/

#include <ctime>
#include <cmath>

#include <itpp/base/itassert.h>
#include <itpp/base/vec.h>
#include <itpp/base/scalfunc.h>

using std::cout;
using std::endl;
using std::abs;
using itpp::it_error_f; // needed for macros it_error_if and it_error

#ifdef _MSC_VER
double lgamma(double x)
{
	it_error_if(x<0,"lgamma: only defined for x>=0");
	return (x+0.5)*log(x+5.5)-x-5.5+log((2.50662827510730+190.9551718930764/(x+1)-216.8366818437280/(x+2)+60.19441764023333/(x+3)-3.08751323928546/(x+4)+0.00302963870525/(x+5)-0.00001352385959072596/(x+6))/x);
}

namespace itpp { 
double gamma(double x)
{
	double s=(2.50662827510730+190.9551718930764/(x+1)-216.8366818437280/(x+2)+60.19441764023333/(x+3)-3.08751323928546/(x+4)+0.00302963870525/(x+5)-0.00001352385959072596/(x+6))/x;
	if (s<0) return -exp((x+0.5)*log(x+5.5)-x-5.5+log(-s));
	else return exp((x+0.5)*log(x+5.5)-x-5.5+log(s));
}
}

double cbrt(double x)
{
	it_error("cbrt() not yet implemented for MS C++");
    return 0.0;
}

#endif

namespace itpp { 

  // Part of C99.
#ifdef _MSC_VER
  double erfc(double Y)
  {
    int  ISW,I;
    double P[4],Q[3],P1[6],Q1[5],P2[4],Q2[3];
    double XMIN,XLARGE,SQRPI;
    double X,RES,XSQ,XNUM,XDEN,XI,XBIG,ERFCret;
    P[1]=0.3166529;
    P[2]=1.722276;
    P[3]=21.38533;
    Q[1]=7.843746;
    Q[2]=18.95226;
    P1[1]=0.5631696;
    P1[2]=3.031799;
    P1[3]=6.865018;
    P1[4]=7.373888;
    P1[5]=4.318779e-5;
    Q1[1]=5.354217;
    Q1[2]=12.79553;
    Q1[3]=15.18491;
    Q1[4]=7.373961;
    P2[1]=5.168823e-2;
    P2[2]=0.1960690;
    P2[3]=4.257996e-2;
    Q2[1]=0.9214524;
    Q2[2]=0.1509421;
    XMIN=1.0E-5;
    XLARGE=4.1875E0;
    XBIG=9.0;
    SQRPI=0.5641896;
    X=Y;
    ISW=1;
    if (X<0) {
      ISW=-1;
      X=-X;
    }
    if (X<0.477) {
      if (X>=XMIN) {
	XSQ=X*X;
	XNUM=(P[1]*XSQ+P[2])*XSQ+P[3];
	XDEN=(XSQ+Q[1])*XSQ+Q[2];
	RES=X*XNUM/XDEN;
      }
      else RES=X*P[3]/Q[2];
      if (ISW==-1) RES=-RES;
      RES=1.0-RES;
      goto slut;
    }
    if (X>4.0) {
      if (ISW>0) goto ulf;
      if (X<XLARGE) goto eva;
      RES=2.0;
      goto slut;
    }
    XSQ=X*X;
    XNUM=P1[5]*X+P1[1];
    XDEN=X+Q1[1];
    for(I=2;I<=4;I++) {
      XNUM=XNUM*X+P1[I];
      XDEN=XDEN*X+Q1[I];
    }
    RES=XNUM/XDEN;
    goto elin;
  ulf:  	if (X>XBIG) goto fred;
  eva:  	XSQ=X*X;
    XI=1.0/XSQ;
    XNUM=(P2[1]*XI+P2[2])*XI+P2[3];
    XDEN=XI+Q2[1]*XI+Q2[2];
    RES=(SQRPI+XI*XNUM/XDEN)/X;
  elin:	RES=RES*exp(-XSQ);
    if (ISW==-1) RES=2.0-RES;
    goto slut;
  fred:	RES=0.0;
  slut:	ERFCret=RES;
    return  ERFCret;
  }

  double erf(double x)
  {
    return (1.0-erfc(x));
  }
#else

  double gamma(double x)
  {
    double lg = lgamma(x);
    return signgam*exp(lg);

  }
#endif

  double Qfunc(double x)
  {
    return 0.5*erfc(x/1.41421356237310);
  }

  double erfinv(double P)
  {
    double	Y,A,B,X,Z,W,WI,SN,SD,F,Z2,SIGMA;
    double	A1=-.5751703,A2=-1.896513,A3=-.5496261E-1;
    double	B0=-.1137730,B1=-3.293474,B2=-2.374996,B3=-1.187515;
    double	C0=-.1146666,C1=-.1314774,C2=-.2368201,C3=.5073975e-1;
    double	D0=-44.27977,D1=21.98546,D2=-7.586103;
    double	E0=-.5668422E-1,E1=.3937021,E2=-.3166501,E3=.6208963E-1;
    double	F0=-6.266786,F1=4.666263,F2=-2.962883;
    double	G0=.1851159E-3,G1=-.2028152E-2,G2=-.1498384,G3=.1078639E-1;
    double	H0=.9952975E-1,H1=.5211733,H2=-.6888301E-1;
    //	double	RINFM=1.7014E+38;

    X=P;
    SIGMA=sgn(X);
    it_error_if(X<-1 || X>1,"erfinv : argument out of bounds");
    Z=fabs(X);
    if (Z>.85) {
      A=1-Z;
      B=Z;
      W=sqrt(-log(A+A*B));
      if (W>=2.5) {
	if (W>=4.) {
	  WI=1./W;
	  SN=((G3*WI+G2)*WI+G1)*WI;
	  SD=((WI+H2)*WI+H1)*WI+H0;
	  F=W+W*(G0+SN/SD);
	} else {
	  SN=((E3*W+E2)*W+E1)*W;
	  SD=((W+F2)*W+F1)*W+F0;
	  F=W+W*(E0+SN/SD);
	}
      } else {
	SN=((C3*W+C2)*W+C1)*W;
	SD=((W+D2)*W+D1)*W+D0;
	F=W+W*(C0+SN/SD);
      }
    } else {
      Z2=Z*Z;
      F=Z+Z*(B0+A1*Z2/(B1+Z2+A2/(B2+Z2+A3/(B3+Z2))));
    }
    Y=SIGMA*F;
    return Y;
  }


#if !defined(__GLIBC__) || __GLIBC__ < 2
  double asinh(double x)
  {
    return ((x>=0) ? log(x+sqrt(x*x+1)):-log(-x+sqrt(x*x+1)));
  }

  double acosh(double x)
  {
    it_error_if(x<1,"acosh(x): x<1");
    return log(x+sqrt(x*x-1));
  }

  double atanh(double x)
  {
    it_error_if(fabs(x)>=1,"atanh(x): abs(x)>=1");
    return 0.5*log((x+1)/(x-1));
  }

#endif

  //Calculates factorial coefficient for index <= 170.
  double fact(int index)
  {
    it_error_if(index > 170,"\nThe function double factfp(int index) overflows if index > 170. \nUse your head instead!");
    it_error_if(index <  0,"\nThe function double factfp(int index) cannot evaluate if index. < 0");
    double prod = 1;
    for (int i=1; i<=index; i++)
      prod *= double(i);
    return prod;
  }

  long mod(long k, long n)
  {
    if (n==0) {
      return k;
    } else {
      return (k - n * long(floor(double(k)/double(n))) );
    }
  }

  long gcd(long a, long b)
  {
    long v, u, t, q;

    it_assert(a>=0,"long gcd(long a, long b): a and b must be non-negative integers");
    it_assert(b>=0,"long gcd(long a, long b): a and b must be non-negative integers");

    u = std::abs(a);
    v = std::abs(b);
    while (v>0) {
      q = u / v;
      t = u - v*q;
      u = v;
      v = t;
    }
    return(u);
  }


  // Calculates binomial coefficient "n over k".
  double binom(int n, int k) {
    it_error_if(k>n,"Error in double binom(int n, int k).\nn must be larger than k.");
    k = n-k<k ? n-k : k;

    vec talj(k), namn(k);
    int i;
    double out = 1.0;
    for (i=0; i<k; i++)
      namn(i) = double(i+1);

    int pos = 0;
    for (i=n; i>=(n-k+1); i--) {
      talj(pos) = double(i);
      pos++;
    }

    for (i=0; i<k; i++) {
      out *= talj(i) / namn(k-1-i);
    }
    return ( out );
  }

  int binom_i(int n, int k)
  {
    ivec v(n);
    int i, j;

    if (n > (k+1)/2)
      n = k+1-n;

    v = 0;
    v(0) = 1;

    for (i=0; i<k-n; i++)
      for (j=1; j<n; j++)
	v(j) += v(j-1);

    return v(n-1);
  }

  // Calculates the base 10-logarithm of the binomial coefficient "n over k".
  double log_binom(int n, int k) {
    it_error_if(k>n,"Error in double log_binom(int n, int k).\nn must be larger than k.");
    k = n-k<k ? n-k : k;

    int i;
    double out = 0.0;
    for (i=0; i<k; i++)
      out += log10((double)(n-i)) - log10((double)(i+1));

    return out;
  }

} //namespace itpp
