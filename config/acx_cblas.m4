dnl @synopsis ACX_CBLAS
dnl @author Adam Piatyszek <ediap@users.sourceforge.net>
dnl @version 2007-02-15
dnl
dnl This macro looks for a library that implements the CBLAS
dnl interface for BLAS (see http://www.netlib.org/blas/). On
dnl success, it sets the CBLAS_LIBS output variable to hold the
dnl requisite library linkages. Besides, it defines HAVE_CBLAS.
dnl
dnl To link with CBLAS, you should link with:
dnl
dnl     $CBLAS_LIBS $BLAS_LIBS $LIBS
dnl
dnl in that order. BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.
dnl
dnl The user may also use --with-cblas=<lib> in order to use some
dnl specific CBLAS library <lib>. In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as was
dnl used to compile the BLAS library.
dnl
dnl @license GPLWithACException

AC_DEFUN([ACX_CBLAS], [
AC_REQUIRE([ACX_BLAS])

AC_ARG_WITH(cblas,
  [AC_HELP_STRING([--with-cblas=<lib>], [use CBLAS library <lib>])])
case $with_cblas in
  yes | "") ;;
  no) acx_cblas_ok=disabled ;;
  -* | */* | *.a | *.so | *.so.* | *.o) CBLAS_LIBS="$with_cblas" ;;
  *) CBLAS_LIBS="-l$with_cblas" ;;
esac

test "$acx_cblas_ok" != disabled && acx_cblas_ok=no
cblas_acml_ok=no

# We cannot use CBLAS if BLAS is not found
test "$acx_blas_ok" != yes && acx_cblas_ok=noblas

# First, check CBLAS_LIBS environment variable
if test "x$CBLAS_LIBS" != x; then
  save_LIBS="$LIBS"; LIBS="$CBLAS_LIBS $BLAS_LIBS $LIBS"
  AC_MSG_CHECKING([for cblas_sgemm in $CBLAS_LIBS])
  AC_TRY_LINK_FUNC(cblas_sgemm, [acx_cblas_ok=yes])
  AC_MSG_RESULT($acx_cblas_ok)
  # Special check for ACML CBLAS
  if test "$acx_cblas_ok" = no; then
    AC_MSG_CHECKING([for sgemm in $CBLAS_LIBS])
    AC_TRY_LINK_FUNC(sgemm, [acx_cblas_ok=yes], [CBLAS_LIBS=""])
    AC_MSG_RESULT($acx_cblas_ok)
    if test "$acx_cblas_ok" = yes; then
      AC_CHECK_HEADER([acml.h], [cblas_acml_ok=yes], [acx_cblas_ok=no])
    fi
  fi
  LIBS="$save_LIBS"
fi

# CBLAS linked to by default?  (it is sometimes included in BLAS)
if test "$acx_cblas_ok" = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS"
  AC_CHECK_FUNC(cblas_sgemm, [acx_cblas_ok=yes])
  LIBS="$save_LIBS"
fi

# CBLAS from ACML linked to by default using FLIBS?
if test "$acx_cblas_ok" = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS"
  AC_CHECK_FUNC(sgemm, [acx_cblas_ok=yes])
  if test "$acx_cblas_ok" = yes; then
    AC_CHECK_HEADER([acml.h], [cblas_acml_ok=yes], [acx_cblas_ok=no])
  fi
  LIBS="$save_LIBS"
fi

# Generic CBLAS library?
for cblas in cblas gslcblas; do
  if test "$acx_cblas_ok" = no; then
    save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
    AC_CHECK_LIB($cblas, cblas_sgemm, [acx_cblas_ok=yes; CBLAS_LIBS="-l$cblas"])
    LIBS="$save_LIBS"
  fi
done

AC_SUBST(CBLAS_LIBS)

# Finally, define HAVE_CBLAS
if test "$acx_cblas_ok" = yes; then
  AC_DEFINE(HAVE_CBLAS, 1, [Define if you have CBLAS library.])
  if test "$cblas_acml_ok" = yes; then
    AC_DEFINE(HAVE_CBLAS_ACML, 1, [Define if you have ACML CBLAS library.])
  fi
fi

])dnl ACX_CBLAS
