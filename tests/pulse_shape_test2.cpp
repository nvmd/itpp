/*!
 * \file 
 * \brief Pulse shaping classes test program
 * \author Adam Piatyszek
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

#include <iomanip>
#include <itpp/itcomm.h>

using namespace std;
using namespace itpp;

int main(int argc, char **argv)
{
  // set precision
  cout.setf(ios::scientific);
  cout.precision(16);

  Parser p;
  QPSK qpsk;
  Root_Raised_Cosine<complex<double> > rrc_tx(0.5), rrc_rx(0.5);
  Raised_Cosine<complex<double> > rc_tx(0.5);

  vec ref_rrc_pulse, ref_rc_pulse;
  cvec ref_symbols, ref_rrc_samples, ref_rc_samples, ref_rec_symbols;
  int error_status = 0;

  // generate tested data
  vec rrc_pulse = rrc_tx.get_pulse_shape();
  vec rc_pulse = rc_tx.get_pulse_shape();
  cvec symbols = qpsk.modulate_bits(randb(100));
  cvec rrc_samples = rrc_tx.shape_symbols(symbols);
  cvec rec_symbols = rrc_rx.shape_samples(rrc_samples);
  cvec rc_samples = rc_tx.shape_symbols(symbols);

  if (argc == 1) {
    // print reference data to stdout
    cout << "rrc_pulse = " << rrc_pulse << endl;
    cout << "rc_pulse = " << rc_pulse << endl;
    cout << "symbols = " << symbols << endl;
    cout << "rrc_samples = " << rrc_samples << endl;
    cout << "rc_samples = " << rc_samples << endl;
    cout << "rec_symbols = " << rec_symbols << endl;
  } 
  else {
    // parse reference data
    p.set_silentmode();
    p.init(argv[1]);
    p.get(ref_rrc_pulse, "rrc_pulse");
    p.get(ref_rc_pulse, "rc_pulse"); 
    p.get(ref_symbols, "symbols");
    p.get(ref_rrc_samples, "rrc_samples");
    p.get(ref_rc_samples, "rc_samples");
    p.get(ref_rec_symbols, "rec_symbols");
    // compare tested data with reference data
    error_status += static_cast<int>(!within_tolerance(rrc_pulse, 
						       ref_rrc_pulse));
    error_status += static_cast<int>(!within_tolerance(rc_pulse, ref_rc_pulse));
    error_status += static_cast<int>(!within_tolerance(symbols, ref_symbols));
    error_status += static_cast<int>(!within_tolerance(rrc_samples, 
						       ref_rrc_samples));
    error_status += static_cast<int>(!within_tolerance(rc_samples,
						       ref_rc_samples));
    error_status += static_cast<int>(!within_tolerance(rec_symbols,
						       ref_rec_symbols));
    if (error_status > 0)
      cerr << "*** pulse_shape_test error status = " << error_status 
	   << " ***" << endl;
  }
  return error_status;
}
