dnl @synopsis ACX_CBLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the CBLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/). On
dnl success, it sets the CBLAS_LIBS output variable to hold the
dnl requisite library linkages.
dnl
dnl To link with CBLAS, you should link with:
dnl
dnl 	$CBLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_CBLAS), and
dnl is sometimes necessary in order to link with F77 libraries. Users
dnl will also need to use AC_F77_DUMMY_MAIN (see the autoconf manual),
dnl for the same reason.
dnl
dnl ATLAS libraries are also searched for. The user may also use 
dnl --with-blas=<lib> in order to use some specific CBLAS library <lib>. 
dnl In order to link successfully, however, be aware that you will 
dnl probably need to use the same Fortran compiler (which can be set 
dnl via the F77 env. var.) as was used to compile the CBLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a CBLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands to
dnl run it if it is not found. If ACTION-IF-FOUND is not specified, the
dnl default action will define HAVE_CBLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @category InstalledPackages
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @version 2001-12-13
dnl @license GPLWithACException
dnl @modified 2005-07-04 by Adam Piatyszek <ediap@users.sourceforge.net>

AC_DEFUN([ACX_CBLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_cblas_ok=no

AC_ARG_WITH(cblas,
	[AC_HELP_STRING([--with-cblas=<lib>], [use CBLAS library <lib>])])
case $with_cblas in
	yes | "") ;;
	no) acx_cblas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) CBLAS_LIBS="$with_cblas" ;;
	*) CBLAS_LIBS="-l$with_cblas" ;;
esac

# Get fortran linker names of CBLAS functions to check for.
AC_F77_FUNC(sgemm)

acx_cblas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check CBLAS_LIBS environment variable
if test $acx_cblas_ok = no; then
if test "x$CBLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$CBLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for cblas_sgemm in $CBLAS_LIBS])
	AC_TRY_LINK_FUNC(cblas_sgemm, [acx_cblas_ok=yes], [CBLAS_LIBS=""])
	AC_MSG_RESULT($acx_cblas_ok)
	LIBS="$save_LIBS"
fi
fi

# CBLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_cblas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[acx_cblas_ok=yes
			 CBLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# CBLAS-ATLAS library in Gentoo?
if test $acx_cblas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[acx_cblas_ok=yes
			 CBLAS_LIBS="-lcblas -lblas -latlas"],
			[], [-lblas -latlas])],
			[], [-latlas])])
fi

# Generic CBLAS library?
if test $acx_cblas_ok = no; then
	AC_CHECK_LIB(cblas, cblas_sgemm, 
		[acx_cblas_ok=yes; CBLAS_LIBS="-lcblas -lblas"], 
		[], [-lblas])
fi

AC_SUBST(CBLAS_LIBS)

LIBS="$acx_cblas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_cblas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_CBLAS,1,[Define if you have a CBLAS library.]),[$1])
        :
else
        acx_cblas_ok=no
        $2
fi
])dnl ACX_CBLAS
