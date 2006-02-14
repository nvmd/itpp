dnl @synopsis ACX_BLAS
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @author Adam Piatyszek <ediap@users.sourceforge.net>
dnl @version 2006-01-31
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/). On
dnl success, it sets the BLAS_LIBS output variable to hold the
dnl requisite library linkages. Besides, it defines HAVE_BLAS.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS), and
dnl is sometimes necessary in order to link with F77 libraries. Users
dnl will also need to use AC_F77_DUMMY_MAIN (see the autoconf manual),
dnl for the same reason.
dnl
dnl Many libraries are searched for, e.g. MKL, ACML or ATLAS. The
dnl user may also use --with-blas=<lib> in order to use some specific
dnl BLAS library <lib>. In order to link successfully, however, be
dnl aware that you will probably need to use the same Fortran compiler
dnl (which can be set via the F77 env. var.) as was used to compile the
dnl BLAS library.
dnl
dnl This macro requires autoconf 2.50 or later.

AC_DEFUN([ACX_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])

AC_ARG_WITH(blas,
  [AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
  yes | "") ;;
  no) acx_blas_ok=disabled ;;
  -* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
  *) BLAS_LIBS="-l$with_blas" ;;
esac

test x"$acx_blas_ok" != xdisabled && acx_blas_ok=no

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(sgemm)
AC_F77_FUNC(dgemm)

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
  if test x"$BLAS_LIBS" != x; then
    save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
    AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
    AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
    AC_MSG_RESULT($acx_blas_ok)
    LIBS="$save_LIBS"
  fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS"
  AC_CHECK_FUNC($sgemm, [acx_blas_ok=yes])
  LIBS="$save_LIBS"
fi

# BLAS in MKL library? 
# (http://www.intel.com/cd/software/products/asmo-na/eng/perflib/mkl/index.htm)
blas_mkl_ok=no
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(mkl, $sgemm, 
    [acx_blas_ok=yes; blas_mkl_ok=yes; BLAS_LIBS="-lmkl -lguide -lpthread"], 
      [], [-lguide -lpthread])
fi

# BLAS in ACML library? (http://developer.amd.com/acml.aspx)
blas_acml_ok=no
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(acml, $sgemm, 
    [acx_blas_ok=yes; blas_acml_ok=yes; BLAS_LIBS="-lacml"])
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
blas_atlas_ok=no
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(atlas, ATL_xerbla,	
    [for blas in f77blas blas; do 
      AC_CHECK_LIB($blas, $sgemm, 
        [AC_CHECK_LIB(cblas, cblas_dgemm,
          [acx_blas_ok=yes; blas_atlas_ok=yes; 
            BLAS_LIBS="-lcblas -l$blas -latlas"],
          [], [-l$blas -latlas])],
        [], [-latlas])
    done])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(blas, $sgemm,
    [AC_CHECK_LIB(dgemm, $dgemm,
      [AC_CHECK_LIB(sgemm, $sgemm,
        [acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
        [], [-lblas])],
      [], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
  if test "x$GCC" != xyes; then # only works with Sun CC
    AC_CHECK_LIB(sunmath, acosp,
      [AC_CHECK_LIB(sunperf, $sgemm,
        [acx_blas_ok=yes; BLAS_LIBS="-xlic_lib=sunperf -lsunmath"], [], 
          [-lsunmath])])
  fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(complib.sgimath, $sgemm,
    [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(blas, $sgemm,
    [AC_CHECK_LIB(essl, $sgemm,
      [acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"], [], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
  AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"

# Finally, define HAVE_BLAS
if test x"$acx_blas_ok" = xyes; then
  AC_DEFINE(HAVE_BLAS, 1, [Define if you have a BLAS library.])
fi

])dnl ACX_BLAS
