/*!
 * \file
 * \brief Resampling functions - source file
 * \author Tony Ottosson and Adam Piatyszek
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

#include <itpp/signal/resampling.h>


namespace itpp
{

// ----------------------------------------------------------------------
// Instantiations
// ----------------------------------------------------------------------

template vec repeat(const vec &v, int norepeats);
template cvec repeat(const cvec &v, int norepeats);
template svec repeat(const svec &v, int norepeats);
template ivec repeat(const ivec &v, int norepeats);
template bvec repeat(const bvec &v, int norepeats);

template mat repeat(const mat &m, int norepeats);
template cmat repeat(const cmat &m, int norepeats);
template smat repeat(const smat &m, int norepeats);
template imat repeat(const imat &m, int norepeats);
template bmat repeat(const bmat &m, int norepeats);

template vec upsample(const vec &v, int usf);
template cvec upsample(const cvec &v, int usf);
template svec upsample(const svec &v, int usf);
template ivec upsample(const ivec &v, int usf);
template bvec upsample(const bvec &v, int usf);

template mat upsample(const mat &v, int usf);
template cmat upsample(const cmat &v, int usf);
template smat upsample(const smat &v, int usf);
template imat upsample(const imat &v, int usf);
template bmat upsample(const bmat &v, int usf);

template void upsample(const vec &v, int usf,  vec & u);
template void upsample(const cvec &v, int usf,  cvec & u);
template void upsample(const svec &v, int usf,  svec & u);
template void upsample(const ivec &v, int usf,  ivec & u);
template void upsample(const bvec &v, int usf,  bvec & u);

template void upsample(const mat &v, int usf,  mat & u);
template void upsample(const cmat &v, int usf,  cmat & u);
template void upsample(const smat &v, int usf,  smat & u);
template void upsample(const imat &v, int usf,  imat & u);
template void upsample(const bmat &v, int usf,  bmat & u);

template vec lininterp(const vec &v, int usf);
template cvec lininterp(const cvec &v, int usf);

template mat lininterp(const mat &v, int usf);
template cmat lininterp(const cmat &v, int usf);

template void lininterp(const vec &v, int usf,  vec & u);
template void lininterp(const cvec &v, int usf,  cvec & u);

template void lininterp(const mat &v, int usf,  mat & u);
template void lininterp(const cmat &v, int usf,  cmat & u);

template mat lininterp(const mat &m, double f_base, double f_ups, int nrof_samples, double t_start);
template cmat lininterp(const cmat &m, double f_base, double f_ups, int nrof_samples, double t_start);

template vec lininterp(const vec &v, double f_base, double f_ups, int nrof_samples, double t_start);
template cvec lininterp(const cvec &v, double f_base, double f_ups, int nrof_samples, double t_start);

} // namespace itpp
