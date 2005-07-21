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

/*----------------------------------------------------------------------------------------*
 * FastICA for IT++                                                                       *
 *----------------------------------------------------------------------------------------*
 * This code is Copyright (C) 2004 by François CAYRE and Teddy FURON                      *
 *                                    TEMICS Project                                      *
 *                                    INRIA/Rennes (IRISA)                                *
 *                                    Campus Universitaire de Beaulieu                    *
 *                                    35042 RENNES cedex FRANCE                           *
 *                                                                                        *
 * Email : firstname.lastname@irisa.fr                                                    *
 *                                                                                        *
 * This is the IT++ implementation of the original Matlab package FastICA.                * 
 *                                                                                        *
 * Matlab package is Copyright (C) 1998 by Jarmo HURRI, Hugo GÄVERT, Jaakko SÄRELÄ and    *
 *                                         Aapo HYVÄRINEN                                 *
 *                                         Laboratory of Information and Computer Science *
 *                                         Helsinki University of Technology              *
 *                                                                                        *
 * URL : http://www.cis.hut.fi/projects/ica/fastica/about.shtml                           *
 *                                                                                        *
 * If you use results given by this FastICA software in an article for a scientific       *
 * journal, conference proceedings or similar, please include the following original      *
 * reference in the bibliography :                                                        *
 *                                                                                        *
 *     A. Hyvärinen. Fast and Robust Fixed-Point Algorithms for Independent Component     *
 *     Analysis. IEEE Transactions on Neural Networks 10(3):626-634, 1999.                *
 *----------------------------------------------------------------------------------------*
 *                                           DISCLAIMER                                   *
 *                                                                                        *
 * This software package is free software ; you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by the Free Software    * 
 * Foundation ; either version 2 of the License, or any later version.                    *
 *                                                                                        *
 * The software package is distributed in the hope that it will be useful, but WITHOUT    *
 * ANY WARRANTY ; without even the implied warranty of MERCHANTABILITY or FITNESS FOR     * 
 * A PARTICULAR PURPOSE.                                                                  *
 * See the GNU General Public License for more details.                                   *
 *                                                                                        *
 *----------------------------------------------------------------------------------------* 
 * Differences with the original Matlab implementation :                                  * 
 * - no GUI                                                                               *
 * - return something even in the case of a convergence problem                           * 
 * - optimization of SVD decomposition (performed 2 times in Matlab, only 1 time in IT++) *
 * - default approach is SYMM wit non-linearity POW3                                      *
 *----------------------------------------------------------------------------------------*/

#include <itpp/itbase.h>
#include "stdio.h"

using std::cout;
using std::endl;
using namespace itpp;

#ifndef HAVE_LAPACK
#define __THIS_PROGRAM_WILL_NOT_RUN__
#endif

#ifndef HAVE_CBLAS
#define __THIS_PROGRAM_WILL_NOT_RUN__
#endif

#ifndef FASTICA_TEST_FILE
#define __THIS_PROGRAM_WILL_NOT_RUN__
#endif

#ifdef __THIS_PROGRAM_WILL_NOT_RUN__
int main() { 
	cout << "1) LAPACK and CBLAS are needed for this test program" << endl;
	cout << "2) FASTICA_TEST_FILE should be defined" << endl;
}
#else

int main() {
  
  FILE * fpin= NULL; 
  float tmp= 0.0; 

  // Separate nrIC independent components in nrSamples samples
  int nrSamples= 0, nrIC= 0;

    cout << "\nTest program for FastICA / IT++\n"; 
    cout << "\n\nFastICA (C) 1998\nJ. Hurri, H. Gävert, J. Särelä, A. Hyvärinen";
    cout << "\nLaboratory of Information and Computer Science\nHelsinki University of Technology";
    cout << "\n\nIT++ port from Matlab by F. Cayre and T. Furon\nTEMICS Team\nIRISA - INRIA/Rennes" << endl; 

    cout << "\nFastICA is the implementation of : \n"; 
    cout << "A. Hyvärinen, Fast and robust fixed-point algorithms for independent component analysis, IEEE Transactions on Neural Networks 10(3):626-634, 1999.\n" << endl;

    cout << "\nSample test data come from fastiCa C standalone version of fastICA, (C) 2003 Jonathan Marchini : " << endl; 
    cout << "URL : http://www.stats.ox.ac.uk/~marchini/software/fastICA/fastiCa.tgz\n" << endl;
    cout << "==========================================================================" << endl;




  fpin = fopen( FASTICA_TEST_FILE, "r");

  fscanf( fpin, "%d", &nrSamples );
  fscanf( fpin, "%d", &nrIC );

  mat X = zeros( nrIC, nrSamples );

  for ( int i= 0; i< nrSamples; i++ ) 
    for ( int j= 0; j< nrIC; j++ ) {
      fscanf( fpin , "%f", &tmp);
      X(j,i)= tmp; 
    }

  fclose( fpin );

  // Instantiate an ICA object with default parameters : SYMM approach and POW3 non-linearity
  // Be sure that : 
  // - nrSamples = number of samples = nb of columns of the input matrix
  // - nrIC = number of sensors = nb of rows of the input matrix
  Fast_ICA my_fastica( X ); 

  // Set number of independent components to separate : 
  // By default, this value is taken from the dimension of 
  // the input data. This line is for illustration purposes. 
  // May help in some cases. 
  my_fastica.set_nrof_independent_components(nrIC); 

  // Perform ICA
  my_fastica.separate(); 

  cout << "Use default parameters:" << endl;

  // Get results : mixing and separating matrices
  cout << "Mixing matrix = " << my_fastica.get_mixing_matrix() << endl;
  cout << "Separation matrix = " << my_fastica.get_separating_matrix() << endl;

  // Get result : separated independent components
  cout << endl << "separated independent components = " << my_fastica.get_independent_components();

  // Another test with other parameters
  cout << "==========================================================" << endl;
  cout << "Use Gaussian non-linearity and deflation approach :" << endl;

  Fast_ICA my_fastica2(X);

  // Set GAUSS non-linearity
  my_fastica2.set_non_linearity(FICA_NONLIN_GAUSS);

  // Use deflation approach : IC are computed one by one
  my_fastica2.set_approach( FICA_APPROACH_DEFL );

  // Perform ICA
  my_fastica2.separate();

  // Get results
  cout << "Mixing matrix = " << my_fastica2.get_mixing_matrix() << endl;
  cout << "Separation matrix = " << my_fastica2.get_separating_matrix() << endl;
  cout << endl << "separated independent components = " << my_fastica2.get_independent_components() << endl;
  
  cout << "End of Fast_ICA test execution. " << endl; 

  exit( 0 );

}

#endif
