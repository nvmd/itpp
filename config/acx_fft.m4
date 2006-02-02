dnl @synopsis ACX_FFT
dnl @author Adam Piatyszek <ediap@users.sourceforge.net>
dnl @version 2006-01-31
dnl @license GPLWithACException
dnl
dnl This macro looks for some FFT implementation, e.g. FFTW3, MKL, ACML,
dnl etc. On success, it sets the FFT_LIBS output variable to hold the
dnl requisite library linkages. Besides it defines one of the following:
dnl HAVE_FFT, HAVE_FFTW, HAVE_MKL8_FFT, HAVE_ACML_FFT
dnl
dnl The user may also use --with-fft=<lib> in order to use some
dnl specific FFT library <lib>.

AC_DEFUN([ACX_FFT], [

AC_ARG_WITH(fft, [AC_HELP_STRING([--with-fft=<lib>], [use FFT library <lib>])])
case $with_fft in
  yes | "") ;;
  no) acx_fft_ok=disabled ;;
  -* | */* | *.a | *.so | *.so.* | *.o) FFT_LIBS="$with_fft" ;;
  *) FFT_LIBS="-l$with_fft" ;;
esac

test x"$acx_fft_ok" != xdisabled && acx_fft_ok=no

fft_mkl8_ok=no
fft_acml_ok=no
fftw3_ok=no

# First, check FFT_LIBS environment variable
if test "x$FFT_LIBS" != x; then
  save_LIBS="$LIBS"; LIBS="$FFT_LIBS $LIBS"
  AC_MSG_CHECKING([for DftiComputeForward in $FFT_LIBS])
  AC_TRY_LINK_FUNC(DftiComputeForward, 
    [AC_CHECK_HEADER([mkl_dfti.h], [acx_fft_ok=yes; fft_mkl8_ok=yes])])
  AC_MSG_RESULT($fft_mkl8_ok)
  if test $acx_fft_ok = no; then
    AC_MSG_CHECKING([for zfft1dx in $FFT_LIBS])
    AC_TRY_LINK_FUNC(zfft1dx, 
      [AC_CHECK_HEADER([acml.h], [acx_fft_ok=yes; fft_acml_ok=yes], 
        [FFT_LIBS=""])],
      [FFT_LIBS=""])
    AC_MSG_RESULT($fft_acml_ok)
  fi
  if test $acx_fft_ok = no; then
    AC_MSG_CHECKING([for fftw_plan_dft_1d in $FFT_LIBS])
    AC_TRY_LINK_FUNC(fftw_plan_dft_1d, 
      [AC_CHECK_HEADER([fftw3.h], [acx_fft_ok=yes; fftw3_ok=yes], 
        [FFT_LIBS=""])], 
      [FFT_LIBS=""])
    AC_MSG_RESULT($fftw3_ok)
  fi
  LIBS="$save_LIBS"
fi

# FFT in BLAS (MKL) library?
if test "x$acx_fft_ok" = xno; then
  if test "x$blas_mkl_ok" = xyes; then
    save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
    AC_CHECK_FUNC(DftiComputeForward, 
      [AC_CHECK_HEADER([mkl_dfti.h], [acx_fft_ok=yes; fft_mkl8_ok=yes])])
    LIBS="$save_LIBS"
  fi
fi

# FFT in MKL library?
if test "x$acx_fft_ok" = xno; then
  AC_CHECK_LIB(mkl, DftiComputeForward, 
    [AC_CHECK_HEADER([mkl_dfti.h], 
      [acx_fft_ok=yes; fft_mkl8_ok=yes; FFT_LIBS="-lmkl -lguide"])],
    [], [-lguide])
fi

# FFT in BLAS (ACML) library?
if test "x$acx_fft_ok" = xno; then
  if test "x$blas_acml_ok" = xyes; then
    save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
    AC_CHECK_FUNC(zfft1dx, 
      [AC_CHECK_HEADER([acml.h], [acx_fft_ok=yes; fft_acml_ok=yes])])
    LIBS="$save_LIBS"
  fi
fi

# FFT in ACML library?
if test "x$acx_fft_ok" = xno; then
  AC_CHECK_LIB(acml, zfft1dx, 
    [AC_CHECK_HEADER([acml.h], 
      [acx_fft_ok=yes; fft_acml_ok=yes; FFT_LIBS="-lacml"])])
fi

# FFT in FFTW3 library?
if test "x$acx_fft_ok" = xno; then
  AC_CHECK_LIB(fftw3, fftw_plan_dft_1d, 
    [AC_CHECK_HEADER([fftw3.h],
      [acx_fft_ok=yes; fftw3_ok=yes; FFT_LIBS="-lfftw3"])])
fi

# FFT in FFTW3 library (extra -lm)?
if test "x$acx_fft_ok" = xno; then
  AC_CHECK_LIB(fftw3, fftw_plan_dft_1d,
    [AC_CHECK_HEADER([fftw3.h],
      [acx_fft_ok=yes; fftw3_ok=yes; FFT_LIBS="-lfftw3 -lm"])], 
    [], [-lm])
fi

AC_SUBST(FFT_LIBS)

# Finally, define HAVE_*
if test x"$acx_fft_ok" = xyes; then
  AC_DEFINE(HAVE_FFT, 1, [Define if you have FFT library.])
  if test x"$fft_mkl8_ok" = xyes; then
    AC_DEFINE(HAVE_FFT_MKL8, 1, [Define if you have MKL8 FFT library.])
  fi
  if test x"$fft_acml_ok" = xyes; then
    AC_DEFINE(HAVE_FFT_ACML, 1, [Define if you have ACML FFT library.])
  fi
  if test x"$fftw3_ok" = xyes; then
    AC_DEFINE(HAVE_FFTW3, 1, [Define if you have FFTW3 library.])
  fi
fi

])dnl ACX_FFT
