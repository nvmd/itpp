#!/bin/sh
# Run this to generate all the initial makefiles, etc.

DIE=no

check_tool() {
    ($1 --version) </dev/null >/dev/null 2>&1 || {
	echo "*** Error: You must have \"$1\" installed to compile IT++ SVN sources"
        DIE=yes
    }
}

check_tool "autoconf"
check_tool "automake"
check_tool "libtoolize"

test "$DIE" = yes && exit 1

srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.
ORIGDIR=`pwd`
cd "$srcdir" || exit $?


PV="`cat VERSION | cut -d' ' -f1`"
LV="`cat VERSION | cut -d' ' -f2`"
PD="`LC_ALL=C date +\"%B %Y\"`"

sed -e "s/@PACKAGE_VERSION@/${PV}/" -e "s/@LIBRARY_VERSION@/${LV}/" \
	< configure.ac.in > configure.ac || exit;
sed -e "s/@PACKAGE_VERSION@/${PV}/" \
	< itpp.spec.in > itpp.spec || exit;
sed -e "s/@PACKAGE_VERSION@/${PV}/" -e "s/@PACKAGE_DATE@/${PD}/" \
	< itpp-config.1.in > itpp-config.1 || exit;

aclocal -I config || exit;
libtoolize --copy --force --automake || exit;
aclocal -I config || exit;
autoconf || exit;
autoheader || exit;
automake --add-missing --copy || exit;

cd "$ORIGDIR" || exit $?
