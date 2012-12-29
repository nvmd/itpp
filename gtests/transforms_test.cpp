/*!
 * \file
 * \brief Transforms test program
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

#include <itpp/itsignal.h>
#include "gtest/gtest.h"
#include <vector>

#ifdef _OPENMP
//should work on both gcc & msvc (it is an omp standard requirement)
#include <omp.h>
#endif

using namespace itpp;

//set test tolerance (measure of relative and absolute error)
const double max_rel_error = 1e-6;
const double max_abs_error = 1e-6;

//transform results tester
template<typename T>
inline void test_result(const Vec<T>& in, const Vec<T>& ref)
{
  int n = size(in);
  it_assert(n == ref.size(), "compute_rel_error(): input and reference sizes must be equal.");
  for(int i = 0; i < n ; ++i) {
    if(abs(in(i) - ref(i)) < max_abs_error) continue; //handle numbers with absolute value close to zero (relative error can be huge for them)
    double rel_error = abs(in(i) - ref(i)) / abs(in(i));
    ASSERT_LE(rel_error, max_rel_error);
  }
}

//Reference transforms implementations. These function is not intended to be fast.
//They just strictly follows the transform definitions
//reference DFT implementation
template<typename T>
inline cvec ref_dft(const Vec<T>& in)
{
  int n = size(in);
  it_assert(n > 0, "ref_dft(): zero-sized input detected.");
  cvec ret(n);
  for(int i = 0; i < n; ++i) {
    std::complex<double> res = 0.0;
    for(int j = 0; j < n; ++j) {
      res += std::complex<double>(cos(2 * pi * i * j / n), -sin(2 * pi * i * j / n)) * in(j);
    }
    ret(i) = res;
  }
  return ret;
}

//reference IDFT implementation
inline cvec ref_idft(const cvec& in)
{
  int n = size(in);
  it_assert(n > 0, "ref_idft(): zero-sized input detected.");
  cvec ret(n);
  for(int i = 0; i < n; ++i) {
    std::complex<double> res = 0.0;
    for(int j = 0; j < n; ++j) {
      res += std::complex<double>(cos(2 * pi * i * j / n), sin(2 * pi * i * j / n)) * in(j);
    }
    ret(i) = res;
  }
  ret *= 1.0 / n;
  return ret;
}

//Type-II DCT reference implementation
inline vec ref_dct(const vec& in)
{
  int n = size(in);
  it_assert(n > 0, "ref_dct(): zero-sized input detected.");
  vec ret(n);
  for(int i = 0; i < n; ++i) {
    double res = 0.0;
    for(int j = 0; j < n; ++j) {
      res += cos(pi * (j + 0.5) * i / n) * in(j);
    }
    ret(i) = 2 * res;
  }

  // Scale to matlab definition format
  ret /= std::sqrt(2.0 * n);
  ret(0) /= std::sqrt(2.0);

  return ret;
}

//Type-III DCT (IDCT) reference implementation
inline vec ref_idct(const vec& in)
{
  int n = size(in);
  it_assert(n > 0, "ref_dct(): zero-sized input detected.");
  vec tmp = in;
  tmp(0) *= std::sqrt(2.0);
  tmp /= std::sqrt(2.0 * n);

  vec ret(n);
  for(int i = 0; i < n; ++i) {
    double res = 0.0;
    for(int j = 1; j < n; ++j) {
      res += cos(pi * (i + 0.5) * j / n) * tmp(j);
    }
    ret(i) = 2 * res + tmp(0);
  }

  return ret;
}

//----------------------------------------------
//Gtest test cases
//----------------------------------------------

TEST(Transforms, FFTReal)
{
  int N = 16;
  vec x = randn(N);
  cvec y;

  if(!have_fourier_transforms()) FAIL() << "Fourier Transforms are not supported with this library build.";

  //vector processing test
  {
    SCOPED_TRACE("fft_real(x, y) test");
    fft_real(x, y);
    test_result(y, ref_dft(x));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = fft_real(x, N) test, N < length(x)");
    int N_sub = 11; //odd subvector length
    y = fft_real(x, N_sub);
    test_result(y, ref_dft(x(0, N_sub - 1)));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = fft_real(x, N) test, N > length(x)");
    int N_zp = 8; //zero-padding length
    y = fft_real(x, N + N_zp);
    x.set_size(N + N_zp, true);
    test_result(y, ref_dft(x));
  }
}

TEST(Transforms, IFFTReal)
{
  if(!have_fourier_transforms()) FAIL() << "Fourier Transforms are not supported with this library build.";

  //vector processing test
  {
    SCOPED_TRACE("ifft_real(x, y) test");
    int N = 16;
    cvec t = randn_c(N - 1);
    vec y;
    cvec x(N);
    //generate test complex sequence with real spectra
    x.set_subvector(1, 0.5 * (t + conj(reverse(t))));
    x(0) = randn();
    //run transform & test results
    ifft_real(x, y);
    test_result(y, real(ref_idft(x)));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = ifft_real(x, N) test, N < length(x)");
    int N = 16, N_sub = 11; //define odd subvector length to test odd-length transform
    cvec t = randn_c(N - 1);
    vec y;
    cvec x(N);
    //fill subvector samples with Hermitian sequence
    x.set_subvector(1, 0.5 * (t(1, N_sub - 1) + conj(reverse(t(1, N_sub - 1)))));
    x(0) = randn();
    //fill the rest of x with random data (these data should be ignored by IFFT implementation)
    x.set_subvector(N_sub, randn_c(N - N_sub));
    //run transform & test results
    y = ifft_real(x, N_sub);
    test_result(y, real(ref_idft(x(0, N_sub - 1))));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = ifft_real(x, N) test, N > length(x)");
    int N_data = 32, N_zp = 8; //define data and zero-padding length
    cvec t = randn_c(N_data);
    vec y;
    cvec x(N_data + N_zp + 1);
    //generate the test data. sequence posesses Hermitian symmetry after zero-padding with N_zp zeros.
    x(0) = randn();
    x.set_subvector(1, N_zp, std::complex<double>(0));
    x.set_subvector(N_zp + 1, 0.5 * (t + conj(reverse(t))));
    //run transform & test results
    y = ifft_real(x, N_data + 2 * N_zp + 1);
    x.set_size(N_data + 2 * N_zp + 1, true);
    test_result(y, real(ref_idft(x)));
  }
}

TEST(Transforms, FFTCplx)
{
  if(!have_fourier_transforms()) FAIL() << "Fourier Transforms are not supported with this library build.";

  int N = 16;
  cvec x = randn_c(N), y;

  //vector processing test
  {
    SCOPED_TRACE("fft(x, y) test");
    fft(x, y);
    test_result(y, ref_dft(x));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = fft(x, N) test, N < length(x)");
    int N_sub = 11; //odd subvector length
    y = fft(x, N_sub);
    test_result(y, ref_dft(x(0, N_sub - 1)));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = fft(x, N) test, N > length(x)");
    int N_zp = 8; //zero-padding length
    y = fft(x, N + N_zp);
    x.set_size(N + N_zp, true);
    test_result(y, ref_dft(x));
  }
}

TEST(Transforms, IFFTCplx)
{
  if(!have_fourier_transforms()) FAIL() << "Fourier Transforms are not supported with this library build.";

  int N = 16;
  cvec x = randn_c(N), y;

  //vector processing test
  {
    SCOPED_TRACE("ifft(x, y) test");
    ifft(x, y);
    test_result(y, ref_idft(x));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = ifft(x, N) test, N < length(x)");
    int N_sub = 11; //odd subvector length
    y = ifft(x, N_sub);
    test_result(y, ref_idft(x(0, N_sub - 1)));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = ifft(x, N) test, N > length(x)");
    int N_zp = 8; //zero-padding length
    y = ifft(x, N + N_zp);
    x.set_size(N + N_zp, true);
    test_result(y, ref_idft(x));
  }
}

TEST(Transforms, DCT)
{
  if(!have_cosine_transforms()) FAIL() << "Cosine Transforms are not supported with this library build.";

  int N = 16;
  vec x = randn(N), y;

  //vector processing test
  {
    SCOPED_TRACE("dct(x, y) test");
    dct(x, y);
    test_result(y, ref_dct(x));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = dct(x, N) test, N < length(x)");
    int N_sub = 11; //odd subvector length
    y = dct(x, N_sub);
    test_result(y, ref_dct(x(0, N_sub - 1)));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = dct(x, N) test, N > length(x)");
    int N_zp = 8; //zero-padding length
    y = dct(x, N + N_zp);
    x.set_size(N + N_zp, true);
    test_result(y, ref_dct(x));
  }
}

TEST(Transforms, IDCT)
{
  if(!have_cosine_transforms()) FAIL() << "Cosine Transforms are not supported with this library build.";

  int N = 16;
  vec x = randn(N), y;

  //vector processing test
  {
    SCOPED_TRACE("idct(x, y) test");
    idct(x, y);
    test_result(y, ref_idct(x));
  }

  //subvector processing test
  {
    SCOPED_TRACE("y = idct(x, N) test, N < length(x)");
    int N_sub = 11; //odd subvector length
    y = idct(x, N_sub);
    test_result(y, ref_idct(x(0, N_sub - 1)));
  }

  //zero-padded vector processing test
  {
    SCOPED_TRACE("y = idct(x, N) test, N > length(x)");
    int N_zp = 8; //zero-padding length
    y = dct(x, N + N_zp);
    x.set_size(N + N_zp, true);
    test_result(y, ref_dct(x));
  }
}


#ifdef _OPENMP

TEST(Transforms, Multithreading)
{
  if(!have_fourier_transforms()) FAIL() << "Fourier Transforms are not supported with this library build.";

  //set number of threads in the team to the maximum available value
  const int threads_cnt = omp_get_max_threads();
  omp_set_num_threads(threads_cnt);

  //in order to test possible clashes on shared data each thread computes fft of length 128, 256,512.
  //Results are compared with reference implementaion upon exit from the parallel region.
  cvec test_input = randn_c(512);
  std::vector<cvec> outputs_128(threads_cnt);
  std::vector<cvec> outputs_256(threads_cnt);
  std::vector<cvec> outputs_512(threads_cnt);

  #pragma omp parallel
  {
    #pragma omp for
    for(int j = 0; j < threads_cnt; ++j) {
      outputs_128[j] = fft(test_input, 128);
      outputs_256[j] = fft(test_input, 256);
      outputs_512[j] = fft(test_input, 512);
    }
  }
  //check results when single-threaded again
  {
    SCOPED_TRACE("fft 128 results.");

    cvec out_ref_128 = ref_dft(test_input(0, 127));
    for(int j = 0; j < threads_cnt; ++j) {
      test_result(outputs_128[j], out_ref_128);
    }
  }
  {
    SCOPED_TRACE("fft 256 results.");
    cvec out_ref_256 = ref_dft(test_input(0, 255));
    for(int j = 0; j < threads_cnt; ++j) {
      test_result(outputs_256[j], out_ref_256);
    }
  }

  {
    SCOPED_TRACE("fft 512 results.");
    cvec out_ref_512 = ref_dft(test_input);
    for(int j = 0; j < threads_cnt; ++j) {
      test_result(outputs_512[j], out_ref_512);
    }
  }

}
#endif