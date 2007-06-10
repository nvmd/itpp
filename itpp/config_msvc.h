/* itpp/config.h.  Generated from config.h.in by configure.  */
/* itpp/config.h.in.  Generated from configure.ac by autoheader.  */


#ifndef CONFIG_H
#define CONFIG_H


/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if you have the `acosh' function. */
/* #undef HAVE_ACOSH */

/* Define to 1 if you have the `asinh' function. */
/* #undef HAVE_ASINH */

/* Define to 1 if you have the `atanh' function. */
/* #undef HAVE_ATANH */

#if defined(HAVE_ACML) || defined(HAVE_MKL)
/* Define if you have a BLAS library. */
#  define HAVE_BLAS 1
#endif

#if defined(HAVE_ACML) || defined(HAVE_MKL)
/* Define if you have CBLAS library. */
#  define HAVE_CBLAS 1
#endif

#ifdef HAVE_ACML
/* Define if you have ACML CBLAS library. */
#  define HAVE_CBLAS_ACML 1
#endif

/* Define to 1 if you have the `cbrt' function. */
/* #undef HAVE_CBRT */

/* Define to 1 if you have the <cmath> header file. */
#define HAVE_CMATH 1

/* Define to 1 if you have the <complex> header file. */
#define HAVE_COMPLEX 1

/* Define to 1 if you have the declaration of `signgam', and to 0 if you
   don't. */
/* #undef HAVE_DECL_SIGNGAM */

/* Define to 1 if you have the <deque> header file. */
#define HAVE_DEQUE 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the `erf' function. */
/* #undef HAVE_ERF */

/* Define to 1 if you have the `erfc' function. */
/* #undef HAVE_ERFC */

#if defined(HAVE_ACML) || defined(HAVE_MKL)
/* Define if you have FFT library. */
#  define HAVE_FFT 1
#endif

/* Define if you have FFTW3 library. */
/* #undef HAVE_FFTW3 */

#ifdef HAVE_ACML
/* Define if you have ACML FFT library. */
#  define HAVE_FFT_ACML 1
#endif

#ifdef HAVE_MKL
/* Define if you have MKL FFT library. */
#  define HAVE_FFT_MKL 1
#endif

/* Define to 1 if you have the `finite' function. */
/* #undef HAVE_FINITE */

/* Define to 1 if you have the `fpclass' function. */
/* #undef HAVE_FPCLASS */

/* Define to 1 if you have the <fstream> header file. */
#define HAVE_FSTREAM 1

/* Define to 1 if you have the <ieeefp.h> header file. */
/* #undef HAVE_IEEEFP_H */

/* Define to 1 if you have the <inttypes.h> header file. */
/* #undef HAVE_INTTYPES_H */

/* Define to 1 if you have the <iomanip> header file. */
#define HAVE_IOMANIP 1

/* Define to 1 if you have the <iostream> header file. */
#define HAVE_IOSTREAM 1

/* Define to 1 if you have the `isfinite' function. */
/* #undef HAVE_ISFINITE */

/* Define to 1 if you have the `isinf' function. */
/* #undef HAVE_ISINF */

/* Define to 1 if you have the `isnan' function. */
/* #undef HAVE_ISNAN */

#if defined(HAVE_ACML) || defined(HAVE_MKL)
/* Define if you have LAPACK library. */
#  define HAVE_LAPACK 1
#endif

/* Define to 1 if you have the `lgamma' function. */
/* #undef HAVE_LGAMMA */

/* Define to 1 if you have the <limits> header file. */
#define HAVE_LIMITS 1

/* Define to 1 if you have the <list> header file. */
#define HAVE_LIST 1

/* Define to 1 if you have the `log1p' function. */
/* #undef HAVE_LOG1P */

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <queue> header file. */
#define HAVE_QUEUE 1

/* Define to 1 if you have the `rint' function. */
/* #undef HAVE_RINT */

/* Define to 1 if you have the <stdexcept> header file. */
#define HAVE_STDEXCEPT 1

/* Define to 1 if you have the <stdint.h> header file. */
/* #undef HAVE_STDINT_H */

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <string> header file. */
#define HAVE_STRING 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the `tgamma' function. */
/* #undef HAVE_TGAMMA */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define if you want exceptions handling */
/* #undef ITPP_EXCEPTIONS */

/* Name of package */
#define PACKAGE "itpp"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "ediap@users.sourceforge.net"

/* Define to the full name of this package. */
#define PACKAGE_NAME "IT++"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "IT++ 3.99.1"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "itpp"

/* Define to the version of this package. */
#define PACKAGE_VERSION "3.99.1"

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 4

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `unsigned int', as computed by sizeof. */
#define SIZEOF_UNSIGNED_INT 4

/* The size of `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG 4

/* The size of `unsigned long long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG_LONG 8

/* The size of `unsigned short', as computed by sizeof. */
#define SIZEOF_UNSIGNED_SHORT 2

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
/* #undef TIME_WITH_SYS_TIME */

/* Version number of package */
#define VERSION "3.99.1"


#if defined(HAVE_CMATH)
#  include <cmath>
#endif

/* Solaris uses <ieeefp.h> for declaring isnan() and finite() functions */
#if defined(HAVE_IEEEFP_H)
#  include <ieeefp.h>
#endif

/* Microsoft Visual C++ .NET underscore prefixed functions */
#if defined(_MSC_VER)
#  include <cfloat>
#  define HAVE_FINITE 1
#  define finite(x) _finite(x)
#  define HAVE_ISNAN 1
#  define isnan(x) _isnan(x)
#  define HAVE_FPCLASS 1
#  define fpclass(x) _fpclass(x)
#  define HAVE_JN 1
#  define jn(a, b) _jn(a, b)
#  define HAVE_YN 1
#  define yn(a, b) _yn(a, b)
#  define HAVE_J0 1
#  define j0(a) _j0(a)
#  define HAVE_J1 1
#  define j1(a) _j1(a)
#endif /* defined(_MSC_VER) */

#if (! defined(HAVE_ISINF) && defined(HAVE_FPCLASS))
#  define HAVE_ISINF 1
#  define isinf(a) (fpclass(a) == FP_NINF || fpclass(a) == FP_PINF)
#endif

#if (! defined (HAVE_FINITE) && defined (HAVE_ISFINITE))
#  define HAVE_FINITE 1
#  define finite(a) isfinite(a)
#endif

#if (! defined(HAVE_FINITE) && defined(HAVE_ISNAN) && defined(HAVE_ISINF))
#  define HAVE_FINITE 1
#  define finite(a) (! isnan(a) && ! isinf(a))
#endif

#endif /* #ifndef CONFIG_H */

