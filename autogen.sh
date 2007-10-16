#!/bin/sh

[ -f configure.ac.in ] || {
	echo "'./autogen.sh' has to be run in the IT++ source directory";
	echo;
	exit;
}

(libtoolize --version) < /dev/null > /dev/null 2>&1 || {
	echo "You must have libtool installed to compile IT++";
	echo;
	exit;
}

(automake --version) < /dev/null > /dev/null 2>&1 || {
	echo "You must have automake installed to compile IT++";
	echo;
	exit;
}

(autoconf --version) < /dev/null > /dev/null 2>&1 || {
	echo "You must have autoconf installed to compile IT++";
	echo;
	exit;
}

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

