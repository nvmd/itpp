/*!
 * \file 
 * \brief Statistical routines test program
 * \author Tony Ottosson, Adam Piatyszek and Andy Panov
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
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

#include <itpp/itbase.h>
#include <iomanip>

using namespace itpp;
using namespace std;


void display_pdf(Histogram<double>& hist)
{
  cout.setf(ios::fixed);
  
  const int max_asterisks_per_line = 40;
  ivec bins = hist.get_bins();
 
  // compute and display experimental PDF 
  vec exp_pdf = hist.get_pdf();
  double pdf_max = max(exp_pdf);

  cout << "  bin | count |     PDF |" << endl
       << "------+-------+---------+-------------------------------------------"
       << endl;
  for (int i = 0; i < exp_pdf.length(); i++) {
    int num_asterisks = static_cast<int>(exp_pdf(i) * max_asterisks_per_line 
					 / pdf_max);
    cout << setw(5) << setprecision(1) << round_to_zero(hist.get_bin_center(i))
	 << " | " << setw(5) << hist.get_bin(i) << " | "
	 << setw(7) << setprecision(5) << round_to_zero(exp_pdf(i)) << " | ";
    for (int j = 0; j < num_asterisks; j++) {
      cout << "*";
    }
    cout << endl;
  }  
  cout << "------+-------+---------+-------------------------------------------"
       << endl;
  cout << "Histogram trials counter  : " << sum(hist.get_bins()) << endl;
  cout << "Sum of the histogram bins : " << sum(hist.get_bins()) << endl;
  cout << "Sum of PDF values         : " << sum(exp_pdf) << endl;
  cout << "--------------------------------------------------------------------"
       << endl << endl;
}


int main()
{
  cout << "=================================" << endl;
  cout << "  Test of statistical routines   " << endl;
  cout << "=================================" << endl;

  vec a = randn(5);

  cout << "a = " << a << endl << endl;

  cout << "max(a) = " << max(a) << endl;
  cout << "min(a) = " << min(a) << endl;
  cout << "mean(a) = " << mean(a) << endl;
	cout << "geometric_mean(abs(a)) = " << geometric_mean(abs(a)) << endl;
  cout << "norm(a) = " << norm(a) << endl;
  cout << "norm(a, 2) = " << norm(a,2) << endl;
  cout << "norm(a, 1) = " << norm(a,1) << endl;
  cout << "norm(a, \"fro\") = " << norm(a, "fro") << endl;
  cout << "energy(a) = " << energy(a) << endl;
  cout << "variance(a) = " << variance(a) << endl;
  cout << "moment(a, 1) = " << round_to_zero(moment(a, 1)) << endl;
  cout << "moment(a, 2) = " << moment(a, 2) << endl;
  cout << "moment(a, 3) = " << moment(a, 3) << endl;
  cout << "skewness(a) = " << skewness(a) << endl;
  cout << "kurtosisexcess(a) = " << kurtosisexcess(a) << endl;
  cout << "kurtosis(a) = " << kurtosis(a) << endl << endl;

  mat A = randn(5,5);
  cout << "A = " << A << endl << endl;

  cout << "max(A) = " << max(A) << endl;
  cout << "max(A, 1) = " << max(A,1) << endl;
  cout << "max(A, 2) = " << max(A,2) << endl;
  cout << "min(A) = " << min(A) << endl;
  cout << "min(A, 1) = " << min(A,1) << endl;
  cout << "min(A, 2) = " << min(A,2) << endl;
  cout << "mean(A) = " << mean(A) << endl;
	cout << "geometric_mean(abs(A)) = " << geometric_mean(abs(A)) << endl;
  cout << "norm(A) = " << norm(A) << endl;
  cout << "norm(A, 2) = " << norm(A,2) << endl;
  cout << "norm(A, 1) = " << norm(A,1) << endl;
  cout << "norm(A, \"fro\") = " << norm(A, "fro") << endl << endl;

  cout << "=======================" << endl;
  cout << "    Histogram tests    " << endl;
  cout << "=======================" << endl << endl;
 
  // create histogram
  Histogram<double> hist(-3, 3, 21);
 
  // matrix dimension for statistical test
  int mat_dim = 100;
 
  cout << "Experimental PDF of " << mat_dim << "x" << mat_dim
       << " normal distributed random matrix:" << endl << endl;

  // compute histogram for a random matrix
  hist.update(randn(mat_dim, mat_dim));
  display_pdf(hist);

  // reset histogram, so we can start next experiment
  hist.reset(); 

  // compute histogram for a random vector 
  int num_stat_trials = 50000;
 
  cout << "Experimental PDF of "<< num_stat_trials
       << " normal distributed random variables:" << endl << endl;
  
  // compute histogram for random vector
  hist.update(randn(num_stat_trials));
  display_pdf(hist);
 
  // compute CDF. CDF is computed vs. right bin boundaries
  cout << "Experimental CDF (CDF(x) = Pr(a < x), a - random variable) " << endl
       << "of the same vector:" << endl << endl;
  vec exp_cdf = hist.get_cdf();
  for (int i = 0; i < exp_cdf.length(); i++)
    cout << "CDF(" << setw(5) << setprecision(2) << hist.get_bin_right(i) 
	 << ") = " << setw(6) << setprecision(4) << exp_cdf(i) << endl;

  return 0;
}
