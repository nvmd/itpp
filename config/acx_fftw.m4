dnl @synopsis ACX_FFTW([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for FFTW v3 library (see http://www.fftw.org/). On
dnl success, it sets the FFTW_LIBS output variable to hold the
dnl requisite library linkages.
dnl
dnl The user may also use --with-fftw=<lib> in order to use some
dnl specific FFTW library <lib>.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a FFTW
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands to
dnl run it if it is not found. If ACTION-IF-FOUND is not specified, the
dnl default action will define HAVE_FFTW.
dnl
dnl @category InstalledPackages
dnl @author Adam Piatyszek <ediap@users.sourceforge.net>
dnl @version 2006-01-27
dnl @license GPLWithACException

AC_DEFUN([ACX_FFTW], [
acx_fftw_ok=no

AC_ARG_WITH(fftw,
        [AC_HELP_STRING([--with-fftw=<lib>], [use FFTW library <lib>])])
case $with_fftw in
        yes | "") ;;
        no) acx_fftw_ok=disabled ;;
        -* | */* | *.a | *.so | *.so.* | *.o) FFTW_LIBS="$with_fftw" ;;
        *) FFTW_LIBS="-l$with_fftw" ;;
esac

# First, check FFTW_LIBS environment variable
if test "x$FFTW_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$FFTW_LIBS $LIBS"
	AC_MSG_CHECKING([for fftw_plan_dft_1d in $FFTW_LIBS])
        AC_TRY_LINK_FUNC(fftw_plan_dft_1d, [acx_fftw_ok=yes], [FFTW_LIBS=""])
        AC_MSG_RESULT($acx_fftw_ok)
	LIBS="$save_LIBS"
fi

# Generic FFTW library?
if test "x$acx_fftw_ok" = xno; then
	AC_CHECK_LIB(fftw3, fftw_plan_dft_1d, 
		[acx_fftw_ok=yes; FFTW_LIBS="-lfftw3"], [], [-lm])
fi

AC_SUBST(FFTW_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_fftw_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_FFTW,1,[Define if you have FFTW library.]),[$1])
        :
else
#        acx_fftw_ok=no
	:
        $2
fi
])dnl ACX_FFTW
