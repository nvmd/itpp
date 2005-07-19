#include "itpp/itbase.h"

using std::cout;
using std::endl;
using std::complex;
using namespace itpp;

int main(){

  cout << "======================================" << endl;
  cout << "        Test of filter routines" << endl;
  cout << "======================================" << endl;

  // Test signals
  vec x = randn(20), x2 = randn(20);
  cvec cx = randn_c(20), cx2 = randn_c(20);

  cout << "Input signals: " << endl;
  cout << "x = " << x << endl;
  cout << "x2 = " << x2 << endl;
  cout << "cx = " << cx << endl;
  cout << "cx2 = " << cx2 << endl;

  // Filter coefficients
  vec b(10);
  b.ones(); b(0)+= 0.1;
  cvec cb(10);
  cb.ones(); cb(0)+= 0.1;

  cvec ca(2);
  ca(0) = 1.0; ca(1) = -1.0;
  vec a(2);
  a(0) = 1.0; a(1) = -1.0;

  cout << "Filter coefficients: " << endl;
  cout << "b = " << b << endl;
  cout << "cb = " << cb << endl;
  cout << "a = " << a << endl;
  cout << "ca = " << ca << endl;

  vec y, y2;
  vec s1, s2;

  cvec cy, cy2;
  cvec cs1, cs2;
  
  // ---------------------------------------------------------------
  cout << endl << "-----------------------------------------------" << endl;
  cout << "MA Filter: " << endl;

  MA_Filter<double,double,double> H(b);
  MA_Filter<complex<double>,double,complex<double> > CH(b);
  MA_Filter<complex<double>,complex<double>,complex<double> > C(cb);

  y = H(x);
  s1 = H.get_state();
  cout << "y = " << y << endl;
  cout << "s1 = " << s1 << endl;
  y2 = H(x2);
  s2 = H.get_state();
  cout << "y2 = " << y2 << endl;
  cout << "s2 = " << s2 << endl;
  cout << "Redo, reseting to state s1" << endl;
  H.set_state(s1);
  y2 = H(x2);
  s2 = H.get_state();
  cout << "y2 = " << y2 << endl;
  cout << "s2 = " << s2 << endl;


  cy = CH(cx);
  cs1 = CH.get_state();
  cout << "cy = " << cy << endl;
  cout << "cs1 = " << cs1 << endl;
  cy2 = CH(cx2);
  cs2 = CH.get_state();
  cout << "cy2 = " << cy2 << endl;
  cout << "cs2 = " << cs2 << endl;
  cout << "Redo, reseting to state s1" << endl;
  CH.set_state(cs1);
  cy2 = CH(cx2);
  cs2 = CH.get_state();
  cout << "cy2 = " << cy2 << endl;
  cout << "cs2 = " << cs2 << endl;



  cy = C(cx);
  cout << "cy = " << cy << endl;

  y =  filter(b, 1, x);
  cout << "y = " << y << endl;

  cy =  filter(b, 1, cx);
  cout << "cy = " << cy << endl;

  cy =  filter(cb, 1, cx);
  cout << "cy = " << cy << endl;



  // ---------------------------------------------------------------
  cout << endl << "-----------------------------------------------" << endl;
  cout << "AR Filter: " << endl;

  AR_Filter<double,double,double> HAR(a);
  AR_Filter<complex<double>,double,complex<double> > CHAR(a);
  AR_Filter<complex<double>,complex<double>,complex<double> > CAR(ca);

  y = HAR(x);
  s1 = HAR.get_state();
  cout << "y = " << y << endl;
  cout << "s1 = " << s1 << endl;
  y2 = HAR(x2);
  s2 = HAR.get_state();
  cout << "y2 = " << y2 << endl;
  cout << "s2 = " << s2 << endl;
  cout << "Redo, reseting to state s1" << endl;
  HAR.set_state(s1);
  y2 = HAR(x2);
  s2 = HAR.get_state();
  cout << "y2 = " << y2 << endl;
  cout << "s2 = " << s2 << endl;

  cy = CHAR(cx);
  cs1 = CHAR.get_state();
  cout << "cy = " << cy << endl;
  cout << "cs1 = " << cs1 << endl;
  cy2 = CHAR(cx2);
  cs2 = CHAR.get_state();
  cout << "cy2 = " << cy2 << endl;
  cout << "cs2 = " << cs2 << endl;
  cout << "Redo, reseting to state s1" << endl;
  CHAR.set_state(cs1);
  cy2 = CHAR(cx2);
  cs2 = CHAR.get_state();
  cout << "cy2 = " << cy2 << endl;
  cout << "cs2 = " << cs2 << endl;

  cy = CAR(cx);
  cout << "cy = " << cy << endl;

  y =  filter(1,a,x);
  cout << "y = " << y << endl;

  cy =  filter(1, a, cx);
  cout << "cy = " << cy << endl;

  cy =  filter(1, ca, cx);
  cout << "cy = " << cy << endl;



  // ---------------------------------------------------------------
  cout << endl << "-----------------------------------------------" << endl;
  cout << "ARMA Filter: " << endl;

  ARMA_Filter<double,double,double> HARMA(b, a);
  ARMA_Filter<complex<double>,double,complex<double> > CHARMA(b, a);
  ARMA_Filter<complex<double>,complex<double>,complex<double> > CARMA(cb, ca);

  y = HARMA(x);
  s1 = HARMA.get_state();
  cout << "y = " << y << endl;
  cout << "s1 = " << s1 << endl;
  y2 = HARMA(x2);
  s2 = HARMA.get_state();
  cout << "y2 = " << y2 << endl;
  cout << "s2 = " << s2 << endl;
  cout << "Redo, reseting to state s1" << endl;
  HARMA.set_state(s1);
  y2 = HARMA(x2);
  s2 = HARMA.get_state();
  cout << "y2 = " << y2 << endl;
  cout << "s2 = " << s2 << endl;


  cy = CHARMA(cx);
  cs1 = CHARMA.get_state();
  cout << "cy = " << cy << endl;
  cout << "cs1 = " << cs1 << endl;
  cy2 = CHARMA(cx2);
  cs2 = CHARMA.get_state();
  cout << "cy2 = " << cy2 << endl;
  cout << "cs2 = " << cs2 << endl;
  cout << "Redo, reseting to state s1" << endl;
  CHARMA.set_state(cs1);
  cy2 = CHARMA(cx2);
  cs2 = CHARMA.get_state();
  cout << "cy2 = " << cy2 << endl;
  cout << "cs2 = " << cs2 << endl;



  cy = CARMA(cx);
  cout << "cy = " << cy << endl;

  y =  filter(b, a, x);
  cout << "y = " << y << endl;

  cy =  filter(b, a, cx);
  cout << "cy = " << cy << endl;

  cy =  filter(cb, ca, cx);
  cout << "cy = " << cy << endl;

}
