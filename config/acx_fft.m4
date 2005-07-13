## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id$
##

AC_DEFUN([ACX_FFT],
[
acx_fft_ok=no

AC_ARG_WITH(fft, [  --with-fft=ARG    fft support
      ARG=no          disable FFT support
      ARG=fftw3       FFTW3 support through libfftw3 (default)
      ARG=fftw2       FFTW2 support through libfftw],
  [fft=$withval], [fft=fftw3])

case $fft in
  fftw3*)
    fft=3
    if test "x${SINGLE_PRECISION}" != x; then
      fft_func="fftw_plan_dft_1d"
      fft_libs="fftw3f"
    else
      fft_func="fftw_plan_dft_1d"
      fft_libs="fftw3"
    fi
    ;;
  fftw2*)
    fft=2
    if test "x${SINGLE_PRECISION}" != x; then
      AC_MSG_ERROR([Version 2 of FFTW can not be used with single precision])
    else
      fft_func="fftw3d_f77_create_plan"
      fft_libs="fftw dfftw"
    fi
    ;;
  *)
    fft=no
    ;;
esac

if test "${fft}" != "no"; then

  dnl Check if the library was given in the command line
  if test $acx_fft_ok = no; then
    AC_ARG_WITH(fft-lib, [AC_HELP_STRING([--with-fft-lib=<lib>], [use FFT library <lib>])])
    case $with_fft_lib in
      yes | "") ;;
      no) acx_fft_ok=disable ;;
      -* | */* | *.a | *.so | *.so.* | *.o) LIBS_FFT="$with_fft_lib" ;;
      *) LIBS_FFT="-l$with_fft_lib" ;;
    esac
  fi

  acx_fft_save_LIBS="$LIBS"
  LIBS="$LIBS_FFT $LIBS $FLIBS"

  dnl First, check LIBS_FFT environment variable
  if test $acx_fft_ok = no; then
    if test "x$LIBS_FFT" != x; then
      AC_MSG_CHECKING([for $fft_func in $LIBS_FFT])
      AC_TRY_LINK_FUNC($fft_func, [acx_fft_ok=yes], [])
      AC_MSG_RESULT($acx_fft_ok)
    fi
  fi

  dnl FFTW linked to by default?
  if test $acx_fft_ok = no; then
    AC_CHECK_FUNC($fft_func, [acx_fft_ok=yes])
  fi

  dnl search libraries
  for fftl in $fft_libs; do
    if test $acx_fft_ok = no; then
      AC_CHECK_LIB($fftl, $fft_func,
        [acx_fft_ok=yes; LIBS_FFT="$LIBS_FFT -l$fftl"], [], [$FLIBS])
    fi
  done

  dnl if we have fftw2, search also for the real transforms
  if test "${fft}" = "2"; then
    if test $acx_fft_ok = yes; then

      acx_fft_ok=no;

      AC_CHECK_FUNC(rfftw3d_f77_create_plan, [acx_fft_ok=yes])

      for fftl in rfftw drfftw; do
        if test $acx_fft_ok = no; then
          AC_CHECK_LIB($fftl, rfftw3d_f77_create_plan,
            [acx_fft_ok=yes; LIBS_FFT="$LIBS_FFT -l$fftl"], [], [$LIBS_FFT $FLIBS])
        fi
      done
    fi
  fi

fi

AC_SUBST(LIBS_FFT)
LIBS="$acx_fft_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_fft_ok" = xyes; then
  AC_DEFINE_UNQUOTED(HAVE_FFT, [$fft], [FFT library (fftw2 | fftw3)])
  $1
else
  if test $acx_fft_ok != disable; then
    AC_MSG_WARN([Could not find fftw library. 
                *** Will compile without fftw support])
  fi
  $2
fi
])
