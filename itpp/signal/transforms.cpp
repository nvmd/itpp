/*!
 * \file
 * \brief Fourier, Cosine, Hadamard, Walsh-Hadamard, and 2D Hadamard
 *        transforms - source file 
 * \author Tony Ottosson, Thomas Eriksson, Simon Wood and Adam Piatyszek
 * 
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#if defined(HAVE_FFT_MKL8)
#  include <mkl_dfti.h>
#elif defined(HAVE_FFT_ACML)
#  include <acml.h>
#elif defined(HAVE_FFTW3)
#  include <fftw3.h>
#endif

#include <itpp/signal/transforms.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/math/elem_math.h>
#include <itpp/base/converters.h>


namespace itpp { 

#ifdef HAVE_FFT_MKL8

  //---------------------------------------------------------------------------
  // FFT/IFFT based on MKL
  //---------------------------------------------------------------------------

  void fft(const cvec &in, cvec &out) 
  {
    static DFTI_DESCRIPTOR* fft_handle = NULL;
    static int N;

    out.set_size(in.size(), false);
    if (N != in.size()) {
      N = in.size();
      if (fft_handle != NULL) DftiFreeDescriptor(&fft_handle);
      DftiCreateDescriptor(&fft_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, N);
      DftiSetValue(fft_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
      DftiCommitDescriptor(fft_handle);
    }
    DftiComputeForward(fft_handle, (void *)in._data(), out._data());
  }

  void ifft(const cvec &in, cvec &out)
  {
    static DFTI_DESCRIPTOR* fft_handle = NULL;
    static int N;

    out.set_size(in.size(), false);
    if (N != in.size()) {
      N = in.size();
      if (fft_handle != NULL) DftiFreeDescriptor(&fft_handle);
      DftiCreateDescriptor(&fft_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, N);
      DftiSetValue(fft_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
      DftiSetValue(fft_handle, DFTI_BACKWARD_SCALE, 1.0/N);
      DftiCommitDescriptor(fft_handle);
    }
    DftiComputeBackward(fft_handle, (void *)in._data(), out._data());
  }

  void fft_real(const vec &in, cvec &out)
  {
    static DFTI_DESCRIPTOR* fft_handle = NULL;
    static int N;

    out.set_size(in.size(), false);
    if (N != in.size()) {
      N = in.size();
      if (fft_handle != NULL) DftiFreeDescriptor(&fft_handle);
      DftiCreateDescriptor(&fft_handle, DFTI_DOUBLE, DFTI_REAL, 1, N);
      DftiSetValue(fft_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
      DftiCommitDescriptor(fft_handle);
    }
    DftiComputeForward(fft_handle, (void *)in._data(), out._data());

    // Real FFT does not compute the 2nd half of the FFT points because it
    // is redundant to the 1st half. However, we want all of the data so we
    // fill it in. This is consistent with Matlab's functionality
    int istart = ceil_i(in.size() / 2.0);
    int iend = in.size() - 1;
    int idelta = iend - istart + 1;
    out.set_subvector(istart, iend, reverse(conj(out(1, idelta))));
  }

  void ifft_real(const cvec &in, vec &out)
  {
    static DFTI_DESCRIPTOR* fft_handle = NULL;
    static int N;

    out.set_size(in.size(), false);
    if (N != in.size()) {
      N = in.size();
      if (fft_handle != NULL) DftiFreeDescriptor(&fft_handle);
      DftiCreateDescriptor( &fft_handle, DFTI_DOUBLE, DFTI_REAL, 1, N);
      DftiSetValue(fft_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
      DftiSetValue(fft_handle, DFTI_BACKWARD_SCALE, 1.0/N);
      DftiCommitDescriptor(fft_handle);
    }
    DftiComputeBackward(fft_handle, (void *)in._data(), out._data());
  }

  //---------------------------------------------------------------------------
  // DCT/IDCT based on MKL
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
	c(i) *= std::complex<double>(std::cos(pi*i/N/2), std::sin(-pi*i/N/2))
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
	c(i) *= std::complex<double>(std::cos(pi*i/N/2), std::sin(pi*i/N/2))
	  * std::sqrt(2.0 * N);
      }
      for (int i = N - 1; i >= 1; i--) {
	c(c.size() - i) = c(i) * std::complex<double>(std::cos(pi*i/N), 
						      std::sin(-pi*i/N));
      }
      out = ifft_real(c);
      out.set_size(N, true);
    }
  }

#elif defined(HAVE_FFT_ACML)

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
      zfft1dx(0, 1.0, false, N, (doublecomplex *)in._data(), 1, 
	      (doublecomplex *)out._data(), 1, 
	      (doublecomplex *)comm_ptr->_data(), &info);
    }
    zfft1dx(-1, 1.0, false, N, (doublecomplex *)in._data(), 1, 
	    (doublecomplex *)out._data(), 1, 
	    (doublecomplex *)comm_ptr->_data(), &info);
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
      zfft1dx(0, 1.0/N, false, N, (doublecomplex *)in._data(), 1, 
	      (doublecomplex *)out._data(), 1, 
	      (doublecomplex *)comm_ptr->_data(), &info);
    }
    zfft1dx(1, 1.0/N, false, N, (doublecomplex *)in._data(), 1, 
	    (doublecomplex *)out._data(), 1, 
	    (doublecomplex *)comm_ptr->_data(), &info);
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
      dzfft(0, N, out_re._data(), comm_ptr->_data(), &info);
    }
    dzfft(1, N, out_re._data(), comm_ptr->_data(), &info);

    // Normalise output data
    out_re *= factor;

    // Convert the real Hermitian DZFFT's output to the Matlab's complex form
    vec out_im(in.size()); out_im.zeros();
    out.set_size(in.size(), false);
    out_im.set_subvector(1, reverse(out_re(N/2 + 1, N-1)));
    out_im.set_subvector(N/2 + 1, -out_re(N/2 + 1, N-1));
    out_re.set_subvector(N/2 + 1, reverse(out_re(1, (N-1)/2)));
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
    out.set_subvector(0, real(in(0, in.size()/2)));
    out.set_subvector(in.size()/2 + 1, -imag(in(in.size()/2 + 1, in.size()-1)));

    if (N != in.size()) {
      N = in.size();
      factor = 1.0 / std::sqrt(static_cast<double>(N));
      if (comm_ptr != NULL)
	delete comm_ptr;
      comm_ptr = new vec(5 * in.size() + 100);
      zdfft(0, N, out._data(), comm_ptr->_data(), &info);
    }
    zdfft(1, N, out._data(), comm_ptr->_data(), &info);
    out.set_subvector(1, reverse(out(1, N-1)));

    // Normalise output data
    out *= factor;
  }

  //---------------------------------------------------------------------------
  // DCT/IDCT based on ACML
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
	c(i) *= std::complex<double>(std::cos(pi*i/N/2), std::sin(-pi*i/N/2))
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
	c(i) *= std::complex<double>(std::cos(pi*i/N/2), std::sin(pi*i/N/2))
	  * std::sqrt(2.0 * N);
      }
      for (int i = N - 1; i >= 1; i--) {
	c(c.size() - i) = c(i) * std::complex<double>(std::cos(pi*i/N), 
						      std::sin(-pi*i/N));
      }
      out = ifft_real(c);
      out.set_size(N, true);
    }
  }

#elif defined(HAVE_FFTW3)

  //---------------------------------------------------------------------------
  // FFT/IFFT based on FFTW
  //---------------------------------------------------------------------------

  /* Note: The FFTW_UNALIGNED flag is set when calling the FFTW
     library, as memory alignment may change between calls (this was
     done to fix bug 1418707).

     Setting the FFTW_UNALIGNED flag fixes bug 1418707 for the time
     being but it does not result in maximum possible performance in
     all cases. It appears that a solution that achieves maximum
     performance from FFTW in IT++ would require either to (i)
     redefine IT++ memory management functions to use the type of
     memory alignment required by FFTW, (ii) implement a memory
     re-allocation and data copying mechanism in the IT++/FFTW
     interface function to ensure FFTW is called with properly aligned
     data, or (iii) force the IT++/FFTW interface function to recreate
     the FFT plan whenever memory alignment has changed.  None of
     these solutions was found attractive.
     
  */

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
      p = fftw_plan_dft_1d(N, (fftw_complex *)in._data(), (fftw_complex *)out._data(),
			   FFTW_FORWARD, FFTW_ESTIMATE | FFTW_UNALIGNED);
    }

    // compute FFT using the GURU FFTW interface
    fftw_execute_dft(p, (fftw_complex *)in._data(), (fftw_complex *)out._data());
  }

  void ifft(const cvec &in, cvec &out)
  {
    static int N = 0;
    static double inv_N;
    static fftw_plan p = NULL;
    out.set_size(in.size(), false);

    if (N != in.size()) {
      N = in.size();
      inv_N = 1.0/N;
      if (p != NULL)
	fftw_destroy_plan(p); // destroy the previous plan
      // create a new plan
      p = fftw_plan_dft_1d(N, (fftw_complex *)in._data(), (fftw_complex *)out._data(),
			   FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_UNALIGNED);
    }

    // compute IFFT using the GURU FFTW interface
    fftw_execute_dft(p, (fftw_complex *)in._data(), (fftw_complex *)out._data());

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
      if (p!= NULL)
	fftw_destroy_plan(p); //destroy the previous plan

      // create a new plan
      p = fftw_plan_dft_r2c_1d(N, (double *)in._data(), (fftw_complex *)out._data(),
			       FFTW_ESTIMATE | FFTW_UNALIGNED);
    }

    // compute FFT using the GURU FFTW interface
    fftw_execute_dft_r2c(p, (double *)in._data(), (fftw_complex *)out._data());

    // Real FFT does not compute the 2nd half of the FFT points because it
    // is redundant to the 1st half. However, we want all of the data so we
    // fill it in. This is consistent with Matlab's functionality
    int istart = ceil_i(in.size() / 2.0);
    int iend = in.size() - 1;
    int idelta = iend - istart + 1;
    out.set_subvector(istart, iend, reverse(conj(out(1, idelta))));
  }

  void ifft_real(const cvec &in, vec & out)
  {
    static int N = 0;
    static double inv_N;
    static fftw_plan p = NULL;
    out.set_size(in.size(), false);

    if (N != in.size()) {
      N = in.size();
      inv_N = 1.0/N;
      if (p != NULL) 
	fftw_destroy_plan(p); // destroy the previous plan

      // create a new plan
      p = fftw_plan_dft_c2r_1d(N, (fftw_complex *)in._data(), (double *)out._data(), 
			       FFTW_ESTIMATE | FFTW_UNALIGNED);
    }

    // compute IFFT using the GURU FFTW interface
    fftw_execute_dft_c2r(p, (fftw_complex *)in._data(), (double *)out._data());

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
      if (p!= NULL)
	fftw_destroy_plan(p); // destroy the previous plan

      // create a new plan
      p = fftw_plan_r2r_1d(N, (double *)in._data(), (double *)out._data(), 
			   FFTW_REDFT10, FFTW_ESTIMATE | FFTW_UNALIGNED);
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
      p = fftw_plan_r2r_1d(N, (double *)out._data(), (double *)out._data(),
			   FFTW_REDFT01, FFTW_ESTIMATE | FFTW_UNALIGNED);
    }

    // compute FFT using the GURU FFTW interface
    fftw_execute_r2r(p, (double *)out._data(), (double *)out._data());
  }

#else

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

#endif // HAVE_FFT_MKL8 or HAVE_FFT_ACML or HAVE_FFTW3

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


  template <class T>
  Vec<T> dht(const Vec<T> &v)
  {
    Vec<T> ret(v.size());
    dht(v, ret);
    return ret;
  }

  template <class T>
  void bitrv(Vec<T> &out) {
    int i,j,N1,K,N=out.size();
    T TEMP;

    j=0;
    N1=N-1;
    for (i=0;i<N1;i++) {
      if (i<j) {
	TEMP=out[j];
	out[j]=out[i];
	out[i]=TEMP;
      }
      K=N/2;
      while (K<=j) {
	j=j-K;
	K=K/2;
      }
      j=j+K;
    }
  }

  template <class T>
  void dht(const Vec<T> &vin, Vec<T> &vout)
  {
    T t;
    int m,N,l,k,j,ib,i;

    N = vin.size();
    vout.set_size(N);

    m = levels2bits(N);
    it_assert_debug((1<<m)==N, "dht: The vector size must be a power of two!");

    // This step is separated because it copies vin to vout
    for (ib=0; ib<N; ib+=2) {
      vout(ib) = vin(ib) + vin(ib+1);
      vout(ib+1) = vin(ib) - vin(ib+1);
    }
    N /= 2;

    l = 2;    
    for (i=1; i<m; i++) {
      N /= 2;
      ib = 0;
      for (k=0; k<N; k++) {
	for (j=0; j<l; j++) {
	  t = vout(ib+j);
	  vout(ib+j) += vout(ib+j+l);
	  vout(ib+j+l) = t - vout(ib+j+l);
	}
	ib += 2*l;
      }
      l *= 2;
    }

    vout /= T(std::sqrt(double(vin.size())));
  }

  template <class T>
  void self_dht(Vec<T> &v)
  {
    T t;
    int m,N,l=1,k,j,ib,i;

    N = v.size();

    m = levels2bits(N);
    it_assert_debug((1<<m)==N, "dht: The vector size must be a power of two!");

    for (i=0; i<m; i++) {
      N /= 2;
      ib = 0;
      for (k=0; k<N; k++) {
	for (j=0; j<l; j++) {
	  t = v(ib+j);
	  v(ib+j) += v(ib+j+l);
	  v(ib+j+l) = t - v(ib+j+l);
	}
	ib += 2*l;
      }
      l *= 2;
    }

    v /= T(std::sqrt(double(v.size())));
  }

  template <class T>
  Vec<T> dwht(const Vec<T> &v)
  {
    Vec<T> ret(v.size());

    dwht(v, ret);

    return ret;
  }

  template <class T>
  void dwht(const Vec<T> &vin, Vec<T> &vout)
  {
    dht(vin, vout);
    bitrv(vout);
  }


  template <class T>
  void self_dwht(Vec<T> &v)
  {
    self_dht(v);
    bitrv(v);
  }

  template <class T>
  Mat<T> dht2(const Mat<T> &m)
  {
    Mat<T> ret(m.rows(), m.cols());
    Vec<T> v;
    int i;

    for (i=0; i<m.rows(); i++) {
      v = m.get_row(i);
      self_dht(v);
      ret.set_row(i, v);
    }
    for (i=0; i<m.cols(); i++) {
      v = ret.get_col(i);
      self_dht(v);
      ret.set_col(i, v);
    }

    return transpose(ret);
  }

  template <class T>
  Mat<T> dwht2(const Mat<T> &m)
  {
    Mat<T> ret(m.rows(), m.cols());
    Vec<T> v;
    int i;

    for (i=0; i<m.rows(); i++) {
      v = m.get_row(i);
      self_dwht(v);
      ret.set_row(i, v);
    }
    for (i=0; i<m.cols(); i++) {
      v = ret.get_col(i);
      self_dwht(v);
      ret.set_col(i, v);
    }

    return transpose(ret);
  }

  /////////////////////////////
  //
  //  template instantiation
  //
  /////////////////////////////

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
