#!/bin/sh

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

echo n | libtoolize --copy --force || exit;
aclocal || exit;
#autoheader || exit;
automake --add-missing --copy;
autoconf || exit;
automake || exit;
#./configure $@
