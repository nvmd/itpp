/*!
 * \file
 * \brief Fourier, Cosine, Hadamard, Walsh-Hadamard, and 2D Hadamard
 *        transforms - source file
 * \author Tony Ottosson, Thomas Eriksson, Simon Wood and Adam Piatyszek
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

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#if defined(HAVE_FFT_MKL)
namespace mkl
{
#  include <mkl_dfti.h>
#  undef DftiCreateDescriptor
}
#elif defined(HAVE_FFT_ACML)
namespace acml
{
#  include <acml.h>
}
#elif defined(HAVE_FFTW3)
#  include <fftw3.h>
#endif

#include <itpp/signal/transforms.h>

//! \cond

namespace itpp
{

#if defined(HAVE_FFT_MKL)

//---------------------------------------------------------------------------
// FFT/IFFT based on MKL
//---------------------------------------------------------------------------

void fft(const cvec &in, cvec &out)
{
  static mkl::DFTI_DESCRIPTOR* fft_handle = NULL;
  static int N;

  out.set_size(in.size(), false);
  if (N != in.size()) {
    N = in.size();
    if (fft_handle != NULL) mkl::DftiFreeDescriptor(&fft_handle);
    mkl::DftiCreateDescriptor(&fft_handle, mkl::DFTI_DOUBLE, mkl::DFTI_COMPLEX, 1, N);
    mkl::DftiSetValue(fft_handle, mkl::DFTI_PLACEMENT, mkl::DFTI_NOT_INPLACE);
    mkl::DftiCommitDescriptor(fft_handle);
  }
  mkl::DftiComputeForward(fft_handle, (void *)in._data(), out._data());
}

void ifft(const cvec &in, cvec &out)
{
  static mkl::DFTI_DESCRIPTOR* fft_handle = NULL;
  static int N;

  out.set_size(in.size(), false);
  if (N != in.size()) {
    N = in.size();
    if (fft_handle != NULL) mkl::DftiFreeDescriptor(&fft_handle);
    mkl::DftiCreateDescriptor(&fft_handle, mkl::DFTI_DOUBLE, mkl::DFTI_COMPLEX, 1, N);
    mkl::DftiSetValue(fft_handle, mkl::DFTI_PLACEMENT, mkl::DFTI_NOT_INPLACE);
    mkl::DftiSetValue(fft_handle, mkl::DFTI_BACKWARD_SCALE, 1.0 / N);
    mkl::DftiCommitDescriptor(fft_handle);
  }
  mkl::DftiComputeBackward(fft_handle, (void *)in._data(), out._data());
}

void fft_real(const vec &in, cvec &out)
{
  static mkl::DFTI_DESCRIPTOR* fft_handle = NULL;
  static int N;

  out.set_size(in.size(), false);
  if (N != in.size()) {
    N = in.size();
    if (fft_handle != NULL) mkl::DftiFreeDescriptor(&fft_handle);
    mkl::DftiCreateDescriptor(&fft_handle, mkl::DFTI_DOUBLE, mkl::DFTI_REAL, 1, N);
    mkl::DftiSetValue(fft_handle, mkl::DFTI_PLACEMENT, mkl::DFTI_NOT_INPLACE);
    mkl::DftiCommitDescriptor(fft_handle);
  }
  mkl::DftiComputeForward(fft_handle, (void *)in._data(), out._data());

  // Real FFT does not compute the 2nd half of the FFT points because it
  // is redundant to the 1st half. However, we want all of the data so we
  // fill it in. This is consistent with Matlab's functionality
  int istart = ceil_i(in.size() / 2.0);
  int idelta = in.size() - istart;
  out.set_subvector(istart, reverse(conj(out(1, idelta))));
}

void ifft_real(const cvec &in, vec &out)
{
  static mkl::DFTI_DESCRIPTOR* fft_handle = NULL;
  static int N;

  out.set_size(in.size(), false);
  if (N != in.size()) {
    N = in.size();
    if (fft_handle != NULL) mkl::DftiFreeDescriptor(&fft_handle);
    mkl::DftiCreateDescriptor(&fft_handle, mkl::DFTI_DOUBLE, mkl::DFTI_REAL, 1, N);
    mkl::DftiSetValue(fft_handle, mkl::DFTI_PLACEMENT, mkl::DFTI_NOT_INPLACE);
    mkl::DftiSetValue(fft_handle, mkl::DFTI_BACKWARD_SCALE, 1.0 / N);
    mkl::DftiCommitDescriptor(fft_handle);
  }
  mkl::DftiComputeBackward(fft_handle, (void *)in._data(), out._data());
}

#endif // #ifdef HAVE_FFT_MKL


#if defined(HAVE_FFT_ACML)

//---------------------------------------------------------------------------
// FFT/IFFT based on ACML
//---------------------------------------------------------------------------

void fft(const cvec &in, cvec &out)
{
  static int N = 0;
  static cvec *comm_ptr = NULL;
  int info;
  out.set_size(in.size(), false);

  if (N != in.size()) {
    N = in.size();
    if (comm_ptr != NULL)
      delete comm_ptr;
    comm_ptr = new cvec(5 * in.size() + 100);
    acml::zfft1dx(0, 1.0, false, N, (acml::doublecomplex *)in._data(), 1,
                  (acml::doublecomplex *)out._data(), 1,
                  (acml::doublecomplex *)comm_ptr->_data(), &info);
  }
  acml::zfft1dx(-1, 1.0, false, N, (acml::doublecomplex *)in._data(), 1,
                (acml::doublecomplex *)out._data(), 1,
                (acml::doublecomplex *)comm_ptr->_data(), &info);
}

void ifft(const cvec &in, cvec &out)
{
  static int N = 0;
  static cvec *comm_ptr = NULL;
  int info;
  out.set_size(in.size(), false);

  if (N != in.size()) {
    N = in.size();
    if (comm_ptr != NULL)
      delete comm_ptr;
    comm_ptr = new cvec(5 * in.size() + 100);
    acml::zfft1dx(0, 1.0 / N, false, N, (acml::doublecomplex *)in._data(), 1,
                  (acml::doublecomplex *)out._data(), 1,
                  (acml::doublecomplex *)comm_ptr->_data(), &info);
  }
  acml::zfft1dx(1, 1.0 / N, false, N, (acml::doublecomplex *)in._data(), 1,
                (acml::doublecomplex *)out._data(), 1,
                (acml::doublecomplex *)comm_ptr->_data(), &info);
}

void fft_real(const vec &in, cvec &out)
{
  static int N = 0;
  static double factor = 0;
  static vec *comm_ptr = NULL;
  int info;
  vec out_re = in;

  if (N != in.size()) {
    N = in.size();
    factor = std::sqrt(static_cast<double>(N));
    if (comm_ptr != NULL)
      delete comm_ptr;
    comm_ptr = new vec(5 * in.size() + 100);
    acml::dzfft(0, N, out_re._data(), comm_ptr->_data(), &info);
  }
  acml::dzfft(1, N, out_re._data(), comm_ptr->_data(), &info);

  // Normalise output data
  out_re *= factor;

  // Convert the real Hermitian DZFFT's output to the Matlab's complex form
  vec out_im(in.size());
  out_im.zeros();
  out.set_size(in.size(), false);
  out_im.set_subvector(1, reverse(out_re(N / 2 + 1, N - 1)));
  out_im.set_subvector(N / 2 + 1, -out_re(N / 2 + 1, N - 1));
  out_re.set_subvector(N / 2 + 1, reverse(out_re(1, (N - 1) / 2)));
  out = to_cvec(out_re, out_im);
}

void ifft_real(const cvec &in, vec &out)
{
  static int N = 0;
  static double factor = 0;
  static vec *comm_ptr = NULL;
  int info;

  // Convert Matlab's complex input to the real Hermitian form
  out.set_size(in.size());
  out.set_subvector(0, real(in(0, in.size() / 2)));
  out.set_subvector(in.size() / 2 + 1, -imag(in(in.size() / 2 + 1, in.size() - 1)));

  if (N != in.size()) {
    N = in.size();
    factor = 1.0 / std::sqrt(static_cast<double>(N));
    if (comm_ptr != NULL)
      delete comm_ptr;
    comm_ptr = new vec(5 * in.size() + 100);
    acml::zdfft(0, N, out._data(), comm_ptr->_data(), &info);
  }
  acml::zdfft(1, N, out._data(), comm_ptr->_data(), &info);
  out.set_subvector(1, reverse(out(1, N - 1)));

  // Normalise output data
  out *= factor;
}

#endif // defined(HAVE_FFT_ACML)


#if defined(HAVE_FFTW3)

//---------------------------------------------------------------------------
// FFT/IFFT based on FFTW
//---------------------------------------------------------------------------

void fft(const cvec &in, cvec &out)
{
  static int N = 0;
  static fftw_plan p = NULL;
  out.set_size(in.size(), false);

  if (N != in.size()) {
    N = in.size();
    if (p != NULL)
      fftw_destroy_plan(p); // destroy the previous plan
    // create a new plan
    p = fftw_plan_dft_1d(N, (fftw_complex *)in._data(),
                         (fftw_complex *)out._data(),
                         FFTW_FORWARD, FFTW_ESTIMATE);
  }

  // compute FFT using the GURU FFTW interface
  fftw_execute_dft(p, (fftw_complex *)in._data(),
                   (fftw_complex *)out._data());
}

void ifft(const cvec &in, cvec &out)
{
  static int N = 0;
  static double inv_N;
  static fftw_plan p = NULL;
  out.set_size(in.size(), false);

  if (N != in.size()) {
    N = in.size();
    inv_N = 1.0 / N;
    if (p != NULL)
      fftw_destroy_plan(p); // destroy the previous plan
    // create a new plan
    p = fftw_plan_dft_1d(N, (fftw_complex *)in._data(),
                         (fftw_complex *)out._data(),
                         FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  // compute IFFT using the GURU FFTW interface
  fftw_execute_dft(p, (fftw_complex *)in._data(),
                   (fftw_complex *)out._data());

  // scale output
  out *= inv_N;
}

void fft_real(const vec &in, cvec &out)
{
  static int N = 0;
  static fftw_plan p = NULL;
  out.set_size(in.size(), false);

  if (N != in.size()) {
    N = in.size();
    if (p != NULL)
      fftw_destroy_plan(p); //destroy the previous plan

    // create a new plan
    p = fftw_plan_dft_r2c_1d(N, (double *)in._data(),
                             (fftw_complex *)out._data(),
                             FFTW_ESTIMATE);
  }

  // compute FFT using the GURU FFTW interface
  fftw_execute_dft_r2c(p, (double *)in._data(),
                       (fftw_complex *)out._data());

  // Real FFT does not compute the 2nd half of the FFT points because it
  // is redundant to the 1st half. However, we want all of the data so we
  // fill it in. This is consistent with Matlab's functionality
  int offset = ceil_i(in.size() / 2.0);
  int n_elem = in.size() - offset;
  for (int i = 0; i < n_elem; ++i) {
    out(offset + i) = std::conj(out(n_elem - i));
  }
}

void ifft_real(const cvec &in, vec & out)
{
  static int N = 0;
  static double inv_N;
  static fftw_plan p = NULL;
  out.set_size(in.size(), false);

  if (N != in.size()) {
    N = in.size();
    inv_N = 1.0 / N;
    if (p != NULL)
      fftw_destroy_plan(p); // destroy the previous plan

    // create a new plan
    p = fftw_plan_dft_c2r_1d(N, (fftw_complex *)in._data(),
                             (double *)out._data(),
                             FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  }

  // compute IFFT using the GURU FFTW interface
  fftw_execute_dft_c2r(p, (fftw_complex *)in._data(),
                       (double *)out._data());

  out *= inv_N;
}

//---------------------------------------------------------------------------
// DCT/IDCT based on FFTW
//---------------------------------------------------------------------------

void dct(const vec &in, vec &out)
{
  static int N;
  static fftw_plan p = NULL;
  out.set_size(in.size(), false);

  if (N != in.size()) {
    N = in.size();
    if (p != NULL)
      fftw_destroy_plan(p); // destroy the previous plan

    // create a new plan
    p = fftw_plan_r2r_1d(N, (double *)in._data(),
                         (double *)out._data(),
                         FFTW_REDFT10, FFTW_ESTIMATE);
  }

  // compute FFT using the GURU FFTW interface
  fftw_execute_r2r(p, (double *)in._data(), (double *)out._data());

  // Scale to matlab definition format
  out /= std::sqrt(2.0 * N);
  out(0) /= std::sqrt(2.0);
}

// IDCT
void idct(const vec &in, vec &out)
{
  static int N;
  static fftw_plan p = NULL;
  out = in;

  // Rescale to FFTW format
  out(0) *= std::sqrt(2.0);
  out /= std::sqrt(2.0 * in.size());

  if (N != in.size()) {
    N = in.size();
    if (p != NULL)
      fftw_destroy_plan(p); // destroy the previous plan

    // create a new plan
    p = fftw_plan_r2r_1d(N, (double *)out._data(),
                         (double *)out._data(),
                         FFTW_REDFT01, FFTW_ESTIMATE);
  }

  // compute FFT using the GURU FFTW interface
  fftw_execute_r2r(p, (double *)out._data(), (double *)out._data());
}

#endif // defined(HAVE_FFTW3)


#if defined(HAVE_FFT_MKL) || defined(HAVE_FFT_ACML)

//---------------------------------------------------------------------------
// DCT/IDCT based on MKL or ACML
//---------------------------------------------------------------------------

void dct(const vec &in, vec &out)
{
  int N = in.size();
  if (N == 1)
    out = in;
  else {
    cvec c = fft_real(concat(in, reverse(in)));
    c.set_size(N, true);
    for (int i = 0; i < N; i++) {
      c(i) *= std::complex<double>(std::cos(pi * i / N / 2), std::sin(-pi * i / N / 2))
              / std::sqrt(2.0 * N);
    }
    out = real(c);
    out(0) /= std::sqrt(2.0);
  }
}

void idct(const vec &in, vec &out)
{
  int N = in.size();
  if (N == 1)
    out = in;
  else {
    cvec c = to_cvec(in);
    c.set_size(2*N, true);
    c(0) *= std::sqrt(2.0);
    for (int i = 0; i < N; i++) {
      c(i) *= std::complex<double>(std::cos(pi * i / N / 2), std::sin(pi * i / N / 2))
              * std::sqrt(2.0 * N);
    }
    for (int i = N - 1; i >= 1; i--) {
      c(c.size() - i) = c(i) * std::complex<double>(std::cos(pi * i / N),
                        std::sin(-pi * i / N));
    }
    out = ifft_real(c);
    out.set_size(N, true);
  }
}

#endif // defined(HAVE_FFT_MKL) || defined(HAVE_FFT_ACML)


#if !defined(HAVE_FFT)

void fft(const cvec &in, cvec &out)
{
  it_error("FFT library is needed to use fft() function");
}

void ifft(const cvec &in, cvec &out)
{
  it_error("FFT library is needed to use ifft() function");
}

void fft_real(const vec &in, cvec &out)
{
  it_error("FFT library is needed to use fft_real() function");
}

void ifft_real(const cvec &in, vec & out)
{
  it_error("FFT library is needed to use ifft_real() function");
}

void dct(const vec &in, vec &out)
{
  it_error("FFT library is needed to use dct() function");
}

void idct(const vec &in, vec &out)
{
  it_error("FFT library is needed to use idct() function");
}

#endif // !defined(HAVE_FFT)

cvec fft(const cvec &in)
{
  cvec out;
  fft(in, out);
  return out;
}

cvec fft(const cvec &in, const int N)
{
  cvec in2 = in;
  cvec out;
  in2.set_size(N, true);
  fft(in2, out);
  return out;
}

cvec ifft(const cvec &in)
{
  cvec out;
  ifft(in, out);
  return out;
}

cvec ifft(const cvec &in, const int N)
{
  cvec in2 = in;
  cvec out;
  in2.set_size(N, true);
  ifft(in2, out);
  return out;
}

cvec fft_real(const vec& in)
{
  cvec out;
  fft_real(in, out);
  return out;
}

cvec fft_real(const vec& in, const int N)
{
  vec in2 = in;
  cvec out;
  in2.set_size(N, true);
  fft_real(in2, out);
  return out;
}

vec ifft_real(const cvec &in)
{
  vec out;
  ifft_real(in, out);
  return out;
}

vec ifft_real(const cvec &in, const int N)
{
  cvec in2 = in;
  in2.set_size(N, true);
  vec out;
  ifft_real(in2, out);
  return out;
}

vec dct(const vec &in)
{
  vec out;
  dct(in, out);
  return out;
}

vec idct(const vec &in)
{
  vec out;
  idct(in, out);
  return out;
}


// ----------------------------------------------------------------------
// Instantiation
// ----------------------------------------------------------------------

template vec dht(const vec &v);
template cvec dht(const cvec &v);

template void dht(const vec &vin, vec &vout);
template void dht(const cvec &vin, cvec &vout);

template void self_dht(vec &v);
template void self_dht(cvec &v);

template vec dwht(const vec &v);
template cvec dwht(const cvec &v);

template void dwht(const vec &vin, vec &vout);
template void dwht(const cvec &vin, cvec &vout);

template void self_dwht(vec &v);
template void self_dwht(cvec &v);

template mat  dht2(const mat &m);
template cmat dht2(const cmat &m);

template mat  dwht2(const mat &m);
template cmat dwht2(const cmat &m);

} // namespace itpp

//! \endcond
