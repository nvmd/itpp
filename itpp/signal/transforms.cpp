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

//fft library selector
enum LibraryTag {DoNotHaveLibrary, FFTW, ACML, MKL};
#if defined(HAVE_FFT_MKL)

namespace mkl
{
#  include <mkl_dfti.h>
#  include <mkl_service.h>
#  undef DftiCreateDescriptor
}
const LibraryTag ActiveLibrary = MKL;

#elif defined(HAVE_FFT_ACML)

namespace acml
{
#  include <acml.h>
}
const LibraryTag ActiveLibrary = ACML;

#elif defined(HAVE_FFTW3)

#  include <itpp/fftw3.h>
const LibraryTag ActiveLibrary = FFTW;

#else

const LibraryTag ActiveLibrary = DoNotHaveLibrary;

#endif

#include <itpp/signal/transforms.h>

//! \cond

//multithreading mode selector
enum MultithreadingTag {SingleThreaded = 1, OmpThreaded};
template<MultithreadingTag Tag> class MutexImpl;

#ifdef _OPENMP

#include <omp.h>
const MultithreadingTag ThreadingTag = OmpThreaded;
//number of context records kept per transform type
//see comments for TransformPorvider class for futher details
static const int contexts_per_transform_type = 10;
//specialize mutex for multi-threaded code with OMP
template<> class MutexImpl<OmpThreaded>
{
  omp_lock_t _lck;
  //disable copy-construction and assignment
  MutexImpl(const MutexImpl&);
  MutexImpl& operator=(const MutexImpl&);
public:
  MutexImpl() {omp_init_lock(&_lck);}
  ~MutexImpl() {omp_destroy_lock(&_lck);}
  //lock the mutex
  void lock() {omp_set_lock(&_lck);}
  //try to lock. returns true if ownership is taken
  bool try_lock() {return (omp_test_lock(&_lck)) != 0;}
  //unlock
  void unlock() {omp_unset_lock(&_lck);}
};


#else

const MultithreadingTag ThreadingTag = SingleThreaded;
//number of context records kept per transform type
//see comments for TransformPorvider class for futher details
static const int contexts_per_transform_type = 1;

//specialize mutex for single-threaded code
template<>
class MutexImpl<SingleThreaded>
{
  //disable copy-construction and assignment
  MutexImpl(const MutexImpl&);
  MutexImpl& operator=(const MutexImpl&);
public:
  MutexImpl() {}
  ~MutexImpl() {}
  void lock() {}
  bool try_lock() {return true;}
  void unlock() {}
};

#endif


//define library-specific mutex implementation
typedef MutexImpl<ThreadingTag> Mutex;

//mutex-based lock
class Lock
{
  Mutex& _m;
  //disable copy-construction and assignment
  Lock(const Lock&);
  Lock& operator=(const Lock&);
public:
  Lock(Mutex& m): _m(m) {_m.lock();}
  ~Lock() {_m.unlock();}
};



namespace itpp
{

//define transform types: FFT complex, FFT real, IFFT Complex, IFFT Real, DCT, IDCT
enum TransformType {FFTCplx = 1, FFTReal, IFFTCplx, IFFTReal, DCT, IDCT};
//generic transforms implementation based on transofrm type and specific FFT library
template<TransformType T, LibraryTag L> class Transform;

//FFT library initializer & deinitializer based on mutithreading model
template<MultithreadingTag M, LibraryTag L> void init_fft_library();

#if defined(HAVE_FFT_MKL)
//MKL-specific implementations

//MKL FFT-related notes:
//If multithreading is enabled on ITPP level (and in user's code) MKL FFT descriptors can be created, committed and freed by multiple threads
//In this case Intel recommends to use single-threaded FFT implementation:
//1. "Intel® MKL 10.x threading", http://software.intel.com/en-us/articles/intel-mkl-10x-threading,
//2. "Examples of Using Multi-Threading for FFT Computation", http://software.intel.com/sites/products/documentation/hpc/mkl/mklman/GUID-00422EBE-93C3-4BC9-A621-9BF0A0E93888.htm
//Based on examples, provided by Intel, it seems to be safe to create/commit/free and run FFT on per-thread descriptor
//without additional locking

template<> inline void init_fft_library<SingleThreaded, MKL>() {} //assume no actions required. ITPP does not use threading, so FFT library is free to  use it's own threading implementation
template<> inline void init_fft_library<OmpThreaded, MKL>()
{
  //switch FFT domain of MKL to single-threaded mode as Intel suggests
  //this should work starting from MKL 10.0
  mkl::mkl_domain_set_num_threads(1, MKL_FFT);
}

//---------------------------------------------------------------------------
// FFT/IFFT based on MKL
//---------------------------------------------------------------------------

inline void release_descriptor(mkl::DFTI_DESCRIPTOR* h)
{
  if(h != NULL) {
    if(!mkl::DftiFreeDescriptor(&h))
      it_error("MKL library release_descriptor() failed on DftiFreeDescriptor.");
  }
}

template<> class Transform<FFTCplx, MKL>
{
  mkl::DFTI_DESCRIPTOR* _h;
  int _transform_length;
public:
  Transform(): _h(NULL), _transform_length(0) {}
  typedef cvec InType;
  typedef cvec OutType;
  void compute_transform(const cvec &in, cvec &out) {
    out.set_size(in.size(), false);
    if(_transform_length != in.size()) {
      release_descriptor(_h);
      _transform_length = in.size();

      if(!mkl::DftiCreateDescriptor(&_h, mkl::DFTI_DOUBLE, mkl::DFTI_COMPLEX, 1, _transform_length))
        it_error("MKL library compute_transform() failed on DftiCreateDescriptor.");

      mkl::DftiSetValue(_h, mkl::DFTI_PLACEMENT, mkl::DFTI_NOT_INPLACE);

      if(!mkl::DftiCommitDescriptor(_h))
        it_error("MKL library compute_transform() failed on DftiCommitDescriptor.");
    }
    mkl::DftiComputeForward(_h, (void *)in._data(), out._data());
  }

  void reset() {release_descriptor(_h); *this = Transform();}
};

template<> class Transform<IFFTCplx, MKL>
{
  mkl::DFTI_DESCRIPTOR* _h;
  int _transform_length;
public:
  Transform(): _h(NULL), _transform_length(0) {}
  typedef cvec InType;
  typedef cvec OutType;
  void compute_transform(const cvec &in, cvec &out) {
    out.set_size(in.size(), false);
    if(_transform_length != in.size()) {
      release_descriptor(_h);
      _transform_length = in.size();

      if(!mkl::DftiCreateDescriptor(&_h, mkl::DFTI_DOUBLE, mkl::DFTI_COMPLEX, 1, _transform_length))
        it_error("MKL library compute_transform() failed on DftiCreateDescriptor.");

      mkl::DftiSetValue(_h, mkl::DFTI_PLACEMENT, mkl::DFTI_NOT_INPLACE);
      mkl::DftiSetValue(_h, mkl::DFTI_BACKWARD_SCALE, 1.0 / _transform_length);

      if(!mkl::DftiCommitDescriptor(_h))
        it_error("MKL library compute_transform() failed on DftiCommitDescriptor.");
    }
    mkl::DftiComputeBackward(_h, (void *)in._data(), out._data());
  }

  void reset() {release_descriptor(_h); *this = Transform();}
};

template<> class Transform<FFTReal, MKL>
{
  mkl::DFTI_DESCRIPTOR* _h;
  int _transform_length;
public:
  Transform(): _h(NULL), _transform_length(0) {}
  typedef vec InType;
  typedef cvec OutType;
  void compute_transform(const vec &in, cvec &out) {
    out.set_size(in.size(), false);
    if(_transform_length != in.size()) {
      release_descriptor(_h);
      _transform_length = in.size();

      if(!mkl::DftiCreateDescriptor(&_h, mkl::DFTI_DOUBLE, mkl::DFTI_REAL, 1, _transform_length))
        it_error("MKL library compute_transform() failed on DftiCreateDescriptor.");

      mkl::DftiSetValue(_h, mkl::DFTI_PLACEMENT, mkl::DFTI_NOT_INPLACE);

      if(!mkl::DftiCommitDescriptor(_h))
        it_error("MKL library compute_transform() failed on DftiCommitDescriptor.");
    }
    mkl::DftiComputeForward(_h, (void *)in._data(), out._data());
    // Real FFT does not compute the 2nd half of the FFT points because it
    // is redundant to the 1st half. However, we want all of the data so we
    // fill it in. This is consistent with Matlab's functionality
    int istart = ceil_i(in.size() / 2.0);
    int idelta = in.size() - istart;
    out.set_subvector(istart, reverse(conj(out(1, idelta))));
  }

  void reset() {release_descriptor(_h); *this = Transform();}
};

template<> class Transform<IFFTReal, MKL>
{
  mkl::DFTI_DESCRIPTOR* _h;
  int _transform_length;
public:
  Transform(): _h(NULL), _transform_length(0) {}
  typedef cvec InType;
  typedef vec OutType;
  void compute_transform(const cvec &in, vec &out) {
    out.set_size(in.size(), false);
    if(_transform_length != in.size()) {
      release_descriptor(_h);
      _transform_length = in.size();

      if(!mkl::DftiCreateDescriptor(&_h, mkl::DFTI_DOUBLE, mkl::DFTI_REAL, 1, _transform_length))
        it_error("MKL library compute_transform() failed on DftiCreateDescriptor.");

      mkl::DftiSetValue(_h, mkl::DFTI_PLACEMENT, mkl::DFTI_NOT_INPLACE);
      mkl::DftiSetValue(_h, mkl::DFTI_BACKWARD_SCALE, 1.0 / _transform_length);

      if(!mkl::DftiCommitDescriptor(_h))
        it_error("MKL library compute_transform() failed on DftiCommitDescriptor.");
    }
    mkl::DftiComputeBackward(_h, (void *)in._data(), out._data());
  }

  void reset() {release_descriptor(_h); *this = Transform();}
};

#endif // #ifdef HAVE_FFT_MKL


#if defined(HAVE_FFT_ACML)
//ACML-specific implementations

//ACML FFT-related notes:
//ACML documentation is not very verbose regarding the multithreaded use of the library, but multithreaded ifort-built ACML uses
//OMP internally. AMD recommends linking with SINGLE-THREADED library if OMP is enabled in user's code. Also, they claim that
//single-threaded functions can be used from the multiple threads simultaniously (see http://devgurus.amd.com/thread/141592 for
//multi-threading discussion on AMD dev forums) The thread-safety of library functions is also mentioned in ACML release notes (ver 4.4.0).
//In the following implementation we assume that ACML transform functions can be run simultaneously from different threads safely if they operate
//on different data sets.
template<> inline void init_fft_library<SingleThreaded, ACML>() {} //assume no actions required.
template<> inline void init_fft_library<OmpThreaded, ACML>() {}

//---------------------------------------------------------------------------
// FFT/IFFT based on ACML
//---------------------------------------------------------------------------
template<> class Transform<FFTCplx, ACML>
{
  cvec _scratchpad;
  int _transform_length;
public:
  Transform(): _transform_length(0) {}
  typedef cvec InType;
  typedef cvec OutType;
  void compute_transform(const cvec &in, cvec &out) {
    int info;
    out.set_size(in.size(), false);
    if(_transform_length != in.size()) {
      _transform_length = in.size();

      int min_required_size = 5 * _transform_length + 100; //ACML guides suggest  3*size + 100 here, but ITPP code uses 5.
      if(_scratchpad.size() < min_required_size) _scratchpad.set_size(min_required_size);

      acml::zfft1dx(0, 1.0, false, _transform_length, (acml::doublecomplex *)in._data(), 1,
                    (acml::doublecomplex *)out._data(), 1,
                    (acml::doublecomplex *)_scratchpad._data(), &info);
    }
    acml::zfft1dx(-1, 1.0, false, _transform_length, (acml::doublecomplex *)in._data(), 1,
                  (acml::doublecomplex *)out._data(), 1,
                  (acml::doublecomplex *)_scratchpad._data(), &info);
  }

  void reset() {*this = Transform();}
};

template<> class Transform<IFFTCplx, ACML>
{
  cvec _scratchpad;
  int _transform_length;
public:
  Transform(): _transform_length(0) {}
  typedef cvec InType;
  typedef cvec OutType;
  void compute_transform(const cvec &in, cvec &out) {
    int info;
    out.set_size(in.size(), false);
    if(_transform_length != in.size()) {
      _transform_length = in.size();

      int min_required_size = 5 * _transform_length + 100; //ACML guides suggest  3*size + 100 here, but ITPP code uses 5.
      if(_scratchpad.size() < min_required_size) _scratchpad.set_size(min_required_size);

      acml::zfft1dx(0, 1.0 / _transform_length, false, _transform_length, (acml::doublecomplex *)in._data(), 1,
                    (acml::doublecomplex *)out._data(), 1,
                    (acml::doublecomplex *)_scratchpad._data(), &info);
    }
    acml::zfft1dx(1, 1.0 / _transform_length, false, _transform_length, (acml::doublecomplex *)in._data(), 1,
                  (acml::doublecomplex *)out._data(), 1,
                  (acml::doublecomplex *)_scratchpad._data(), &info);
  }

  void reset() {*this = Transform();}
};

template<> class Transform<FFTReal, ACML>
{
  vec _scratchpad;
  int _transform_length;
public:
  Transform(): _transform_length(0) {}
  typedef vec InType;
  typedef cvec OutType;
  void compute_transform(const vec &in, cvec &out) {
    vec out_re = in;

    int info;
    if(_transform_length != in.size()) {
      _transform_length = in.size();

      int min_required_size = 5 * _transform_length + 100; //ACML guides suggest  3*size + 100 here, but ITPP code uses 5.
      if(_scratchpad.size() < min_required_size) _scratchpad.set_size(min_required_size);

      acml::dzfft(0, _transform_length, out_re._data(), _scratchpad._data(), &info);
    }
    acml::dzfft(1, _transform_length, out_re._data(), _scratchpad._data(), &info);

    // Normalise output data
    double factor = std::sqrt(static_cast<double>(_transform_length));
    out_re *= factor;

    // Convert the real Hermitian DZFFT's output to the Matlab's complex form
    vec out_im(_transform_length);
    out_im(0) = 0.0;
    if(!(_transform_length % 2)) out_im(_transform_length / 2) = 0.0; //even transform length
    out_im.set_subvector(1, reverse(out_re(_transform_length / 2 + 1, _transform_length - 1)));
    out_im.set_subvector(_transform_length / 2 + 1, -out_re(_transform_length / 2 + 1, _transform_length - 1));
    out_re.set_subvector(_transform_length / 2 + 1, reverse(out_re(1, (_transform_length - 1) / 2)));

    out.set_size(_transform_length, false);
    out = to_cvec(out_re, out_im);
  }

  void reset() {*this = Transform();}
};

template<> class Transform<IFFTReal, ACML>
{
  vec _scratchpad;
  int _transform_length;
public:
  Transform(): _transform_length(0) {}
  typedef cvec InType;
  typedef vec OutType;
  void compute_transform(const cvec &in, vec &out) {
    // Convert Matlab's complex input to the real Hermitian form
    out.set_size(in.size());
    out.set_subvector(0, real(in(0, in.size() / 2)));
    out.set_subvector(in.size() / 2 + 1, -imag(in(in.size() / 2 + 1, in.size() - 1)));

    int info;
    if(_transform_length != in.size()) {
      _transform_length = in.size();

      int min_required_size = 5 * _transform_length + 100; //ACML guides suggest  3*size + 100 here, but ITPP code uses 5.
      if(_scratchpad.size() < min_required_size) _scratchpad.set_size(min_required_size);

      acml::zdfft(0, _transform_length, out._data(), _scratchpad._data(), &info);
    }
    acml::zdfft(1, _transform_length, out._data(), _scratchpad._data(), &info);
    out.set_subvector(1, reverse(out(1, _transform_length - 1)));

    // Normalise output data
    double factor = 1.0 / std::sqrt(static_cast<double>(_transform_length));
    out *= factor;
  }

  void reset() {*this = Transform();}
};



#endif // defined(HAVE_FFT_ACML)


#if defined(HAVE_FFTW3)
//FFTW3-specific implementations

//FFTW3-related notes:
//Based on the FFtW3 documentation, it is thread-safe to call fftw_execute family functions simultaniously from several threads assuming that data sets are different in each thread.
//FFTW plans creation-destruction is not thread-safe and should be serialized by the caller. FFTW provides some functions to execute transforms with multiple threads (assuming FFTW
// is compiled and linked with multithreading support). Current ITPP implementation does not use any of them.

template<> inline void init_fft_library<SingleThreaded, FFTW>() {} //assume no actions required.
template<> inline void init_fft_library<OmpThreaded, FFTW>() {}
//define global lock for operations with FFTW plans.
Mutex& get_library_lock()
{
  static Mutex FFTW3LibraryLock;
  return FFTW3LibraryLock;
}
//---------------------------------------------------------------------------
// FFT/IFFT based on FFTW
//---------------------------------------------------------------------------
inline void destroy_plan(fftw_plan p)
{
  if(p != NULL) fftw_destroy_plan(p);  // destroy the plan
}

template<> class Transform<FFTCplx, FFTW>
{
  fftw_plan _p;
  int _transform_length;
public:
  Transform(): _p(NULL), _transform_length(0) {}
  typedef cvec InType;
  typedef cvec OutType;
  void compute_transform(const cvec &in, cvec &out) {
    out.set_size(in.size(), false);
    if(_transform_length != in.size()) {
      Lock l(get_library_lock()); //apply global library lock on plan changes

      _transform_length = in.size();
      destroy_plan(_p); // destroy the previous plan
      // create a new plan (creation of plan guarantees not to return NULL)
      _p = fftw_plan_dft_1d(_transform_length, (fftw_complex *)in._data(),
                            (fftw_complex *)out._data(),
                            FFTW_FORWARD, FFTW_ESTIMATE);
    }
    //compute FFT using the GURU FFTW interface
    fftw_execute_dft(_p, (fftw_complex *)in._data(),
                     (fftw_complex *)out._data());
  }

  void reset() {destroy_plan(_p); *this = Transform();}
};

template<> class Transform<IFFTCplx, FFTW>
{
  fftw_plan _p;
  int _transform_length;
public:
  Transform(): _p(NULL), _transform_length(0) {}
  typedef cvec InType;
  typedef cvec OutType;
  void compute_transform(const cvec &in, cvec &out) {
    out.set_size(in.size(), false);
    if(_transform_length != in.size()) {
      Lock l(get_library_lock()); //apply global library lock on plan changes

      _transform_length = in.size();
      destroy_plan(_p); // destroy the previous plan

      // create a new plan (creation of plan guarantees not to return NULL)
      _p = fftw_plan_dft_1d(_transform_length, (fftw_complex *)in._data(),
                            (fftw_complex *)out._data(),
                            FFTW_BACKWARD, FFTW_ESTIMATE);
    }
    //compute FFT using the GURU FFTW interface
    fftw_execute_dft(_p, (fftw_complex *)in._data(),
                     (fftw_complex *)out._data());
    // scale output
    double inv_N = 1.0 / _transform_length;
    out *= inv_N;

  }

  void reset() {destroy_plan(_p); *this = Transform();}
};

template<> class Transform<FFTReal, FFTW>
{
  fftw_plan _p;
  int _transform_length;
public:
  Transform(): _p(NULL), _transform_length(0) {}
  typedef vec InType;
  typedef cvec OutType;
  void compute_transform(const vec &in, cvec &out) {
    out.set_size(in.size(), false);
    if(_transform_length != in.size()) {
      Lock l(get_library_lock()); //apply global library lock on plan changes

      _transform_length = in.size();
      destroy_plan(_p); // destroy the previous plan

      // create a new plan (creation of plan guarantees not to return NULL)
      _p = fftw_plan_dft_r2c_1d(_transform_length, (double *)in._data(),
                                (fftw_complex *)out._data(),
                                FFTW_ESTIMATE);
    }
    //compute FFT using the GURU FFTW interface
    fftw_execute_dft_r2c(_p, (double *)in._data(),
                         (fftw_complex *)out._data());
    // Real FFT does not compute the 2nd half of the FFT points because it
    // is redundant to the 1st half. However, we want all of the data so we
    // fill it in. This is consistent with Matlab's functionality
    int offset = ceil_i(_transform_length / 2.0);
    int n_elem = _transform_length - offset;
    for(int i = 0; i < n_elem; ++i) {
      out(offset + i) = std::conj(out(n_elem - i));
    }
  }

  void reset() {destroy_plan(_p); *this = Transform();}
};

template<> class Transform<IFFTReal, FFTW>
{
  fftw_plan _p;
  int _transform_length;
public:
  Transform(): _p(NULL), _transform_length(0) {}
  typedef cvec InType;
  typedef vec OutType;
  void compute_transform(const cvec &in, vec &out) {
    out.set_size(in.size(), false);
    if(_transform_length != in.size()) {
      Lock l(get_library_lock()); //apply global library lock on plan changes

      _transform_length = in.size();
      destroy_plan(_p); // destroy the previous plan

      // create a new plan (creation of plan guarantees not to return NULL)
      _p = fftw_plan_dft_c2r_1d(_transform_length, (fftw_complex *)in._data(),
                                (double *)out._data(),
                                FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    }
    //compute FFT using the GURU FFTW interface
    fftw_execute_dft_c2r(_p, (fftw_complex *)in._data(),
                         (double *)out._data());
    // scale output
    double inv_N = 1.0 / _transform_length;
    out *= inv_N;

  }

  void reset() {destroy_plan(_p); *this = Transform();}
};

//---------------------------------------------------------------------------
// DCT/IDCT based on FFTW
//---------------------------------------------------------------------------
template<> class Transform<DCT, FFTW>
{
  fftw_plan _p;
  int _transform_length;
public:
  Transform(): _p(NULL), _transform_length(0) {}
  typedef vec InType;
  typedef vec OutType;
  void compute_transform(const vec &in, vec &out) {
    out.set_size(in.size(), false);
    if(_transform_length != in.size()) {
      Lock l(get_library_lock()); //apply global library lock on plan changes

      _transform_length = in.size();
      destroy_plan(_p); // destroy the previous plan

      // create a new plan (creation of plan guarantees not to return NULL)
      _p = fftw_plan_r2r_1d(_transform_length, (double *)in._data(),
                            (double *)out._data(),
                            FFTW_REDFT10, FFTW_ESTIMATE);
    }
    // compute FFT using the GURU FFTW interface
    fftw_execute_r2r(_p, (double *)in._data(), (double *)out._data());

    // Scale to matlab definition format
    out /= std::sqrt(2.0 * _transform_length);
    out(0) /= std::sqrt(2.0);
  }

  void reset() {destroy_plan(_p); *this = Transform();}
};

template<> class Transform<IDCT, FFTW>
{
  fftw_plan _p;
  int _transform_length;
public:
  Transform(): _p(NULL), _transform_length(0) {}
  typedef vec InType;
  typedef vec OutType;
  void compute_transform(const vec &in, vec &out) {
    out = in;

    // Rescale to FFTW format
    out(0) *= std::sqrt(2.0);
    out /= std::sqrt(2.0 * in.size());

    if(_transform_length != in.size()) {
      Lock l(get_library_lock()); //apply global library lock on plan changes

      _transform_length = in.size();
      destroy_plan(_p); // destroy the previous plan

      // create a new plan (creation of plan guarantees not to return NULL)
      _p = fftw_plan_r2r_1d(_transform_length, (double *)in._data(),
                            (double *)out._data(),
                            FFTW_REDFT01, FFTW_ESTIMATE);
    }
    // compute FFT using the GURU FFTW interface
    fftw_execute_r2r(_p, (double *)out._data(), (double *)out._data());
  }

  void reset() {destroy_plan(_p); *this = Transform();}
};

#endif // defined(HAVE_FFTW3)

#if defined(HAVE_FFT_MKL) || defined(HAVE_FFT_ACML)

//---------------------------------------------------------------------------
// DCT/IDCT based on MKL or ACML
//---------------------------------------------------------------------------

//use FFT on real values to perform DCT
template<> class Transform<DCT, ActiveLibrary>
{
  Transform<FFTReal, ActiveLibrary> _tr;
public:
  Transform() {}
  typedef vec InType;
  typedef vec OutType;
  void compute_transform(const vec &in, vec &out) {
    int N = in.size();
    if(N == 1)
      out = in;
    else {
      cvec c;
      _tr.compute_transform(concat(in, reverse(in)), c);
      c.set_size(N, true);
      for(int i = 0; i < N; i++) {
        c(i) *= std::complex<double>(std::cos(pi * i / N / 2), std::sin(-pi * i / N / 2))
                / std::sqrt(2.0 * N);
      }
      out = real(c);
      out(0) /= std::sqrt(2.0);
    }
  }

  void reset() {_tr.reset();}
};

//use IFFT with real output to perform IDCT
template<> class Transform<IDCT, ActiveLibrary>
{
  Transform<IFFTReal, ActiveLibrary> _tr;
public:
  Transform() {}
  typedef vec InType;
  typedef vec OutType;
  void compute_transform(const vec &in, vec &out) {
    int N = in.size();
    if(N == 1)
      out = in;
    else {
      cvec c = to_cvec(in);
      c.set_size(2 * N, true);
      c(0) *= std::sqrt(2.0);
      for(int i = 0; i < N; i++) {
        c(i) *= std::complex<double>(std::cos(pi * i / N / 2), std::sin(pi * i / N / 2))
                * std::sqrt(2.0 * N);
      }
      for(int i = N - 1; i >= 1; i--) {
        c(c.size() - i) = c(i) * std::complex<double>(std::cos(pi * i / N),
                          std::sin(-pi * i / N));
      }
      _tr.compute_transform(c, out);
      out.set_size(N, true);
    }
  }

  void reset() {_tr.reset();}
};

#endif

#if defined(HAVE_FFT)
//lock-protected transform to serialize accesses to the context from several threads
template<TransformType TT> class LockedTransform : private Transform<TT, ActiveLibrary>
{
  typedef Transform<TT, ActiveLibrary> Base;
  Mutex _m;
public:
  typedef typename Base::InType InType;
  typedef typename Base::OutType OutType;
  LockedTransform() {}
  //release context
  void release_context() {Lock l(_m); Base::reset();}
  void run_transform(const InType& in, OutType& out) {Lock l(_m); Base::compute_transform(in, out);}
};

//Typical multithreaded application creates several threads upon entry to parallel region and join them upon exit from it.
//Threads used to perform parallel computations can either be terminated upon exit from the parallel region or left in the
//parked state, so in the next parallel region application can reuse already created threads from the team instead of the
//time-consuming creation of new threads. There is no way to control threads creation-destruction with OMP and, therefore
//there is no way to implement automatic clean-up of transform computation contexts for each thread (basically, this means that
//we can not appropriately release FFT library resources and this results in memory leak)

//In order solve this problem and implement the FFT transforms in multithreded environment library relyes on the statically
//created pool of transform contexts.

//Each thread willing to run the transfrom queries the context index from transfrom provider. Thread uses assigned index and
//corresponding context to compute the transforms during it's lifetime. Provider assigns contexts in round-robbin fashion, so
//context used by the exited threads are reused by newly created ones.

//Single context can be reused by multiple threads if application created more then contexts_per_transform_type threads performing
//some type of transform.
static bool is_library_initialized = false;

template<TransformType TT> class TransformPorvider
{
  typedef LockedTransform<TT> Transform;
  Transform _transforms[contexts_per_transform_type];
  int _id;
public:
  TransformPorvider(): _id(0) {
    if(!is_library_initialized) {
      //initialize FFT library on first conctruction of any of TransformPorvider objects
      init_fft_library<ThreadingTag, ActiveLibrary>();
      is_library_initialized = true;
    }
  }
  int get_context_id() {
    //assign id in round-robin fashion.
    int ret = _id + 1;
    if(ret == contexts_per_transform_type)
      _id = 0;
    else
      _id = ret;
    return ret;
  }
  void run_transform(int id, const typename Transform::InType& in, typename Transform::OutType& out) {
    _transforms[id - 1].run_transform(in, out);
  }
  //provider destructor. releases context resources.
  //destructor is called after the main() exits, so there is no need to protect context release with mutex
  ~TransformPorvider() {
    for(int i = 0; i < contexts_per_transform_type; ++i)
      _transforms[i].release_context();
  }
};

//TransformPorvider is constructed upon the first request
template<TransformType TT> TransformPorvider<TT>&  get_transform_provider()
{
  static TransformPorvider<TT> p;
  return p;
}

void fft(const cvec &in, cvec &out)
{
  static int context_id = 0;
  #pragma omp threadprivate(context_id)

  if(context_id == 0) {
    //first-time transform call
    #pragma omp critical
    {
      //serialize access to  transform provider to get the id
      context_id = get_transform_provider<FFTCplx>().get_context_id();
    }
  }
  it_assert(in.size() > 0, "fft(): zero-sized input detected");
  //there is no need to serialize here, since provider is constructed at this point
  get_transform_provider<FFTCplx>().run_transform(context_id, in, out);
}

void ifft(const cvec &in, cvec &out)
{
  static int context_id = 0;
  #pragma omp threadprivate(context_id)

  if(context_id == 0) {
    //first-time transform call
    #pragma omp critical
    {
      //serialize access to  transform provider to get the id
      context_id = get_transform_provider<IFFTCplx>().get_context_id();
    }
  }
  it_assert(in.size() > 0, "ifft(): zero-sized input detected");
  //there is no need to serialize here, since provider is constructed at this point
  get_transform_provider<IFFTCplx>().run_transform(context_id, in, out);
}

void fft_real(const vec &in, cvec &out)
{
  static int context_id = 0;
  #pragma omp threadprivate(context_id)

  if(context_id == 0) {
    //first-time transform call
    #pragma omp critical
    {
      //serialize access to  transform provider to get the id
      context_id = get_transform_provider<FFTReal>().get_context_id();
    }
  }
  it_assert(in.size() > 0, "fft_real(): zero-sized input detected");
  //there is no need to serialize here, since provider is constructed at this point
  get_transform_provider<FFTReal>().run_transform(context_id, in, out);
}

void ifft_real(const cvec &in, vec &out)
{
  static int context_id = 0;
  #pragma omp threadprivate(context_id)

  if(context_id == 0) {
    //first-time transform call
    #pragma omp critical
    {
      //serialize access to  transform provider to get the id
      context_id = get_transform_provider<IFFTReal>().get_context_id();
    }
  }
  it_assert(in.size() > 0, "ifft_real(): zero-sized input detected");
  //there is no need to serialize here, since provider is constructed at this point
  get_transform_provider<IFFTReal>().run_transform(context_id, in, out);
}

void dct(const vec &in, vec &out)
{
  static int context_id = 0;
  #pragma omp threadprivate(context_id)

  if(context_id == 0) {
    //first-time transform call
    #pragma omp critical
    {
      //serialize access to  transform provider to get the id
      context_id = get_transform_provider<DCT>().get_context_id();
    }
  }
  it_assert(in.size() > 0, "dct(): zero-sized input detected");
  //there is no need to serialize here, since provider is definitely constructed at this point
  get_transform_provider<DCT>().run_transform(context_id, in, out);
}

void idct(const vec &in, vec &out)
{
  static int context_id = 0;
  #pragma omp threadprivate(context_id)

  if(context_id == 0) {
    //first-time transform call
    #pragma omp critical
    {
      //serialize access to  transform provider to get the id
      context_id = get_transform_provider<IDCT>().get_context_id();
    }
  }
  it_assert(in.size() > 0, "dct(): zero-sized input detected");
  //there is no need to serialize here, since provider is definitely constructed at this point
  get_transform_provider<IDCT>().run_transform(context_id, in, out);
}

bool have_fourier_transforms() {return true;}
bool have_cosine_transforms() {return true;}
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

bool have_fourier_transforms() {return false;}
bool have_cosine_transforms() {return false;}

#endif // defined(HAVE_FFT)


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

vec dct(const vec &in, int N)
{
  vec in2 = in;
  vec out;
  in2.set_size(N, true);
  dct(in2, out);
  return out;
}


vec idct(const vec &in)
{
  vec out;
  idct(in, out);
  return out;
}

vec idct(const vec &in, int N)
{
  vec in2 = in;
  vec out;
  in2.set_size(N, true);
  idct(in2, out);
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
