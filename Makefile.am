include $(top_srcdir)/Makefile.common

ACLOCAL_AMFLAGS = -I config

if HTML_DOC
HTMLDOCDIR = doc
endif

if HAVE_DIFF
TESTDIR = tests
endif

SUBDIRS = itpp $(HTMLDOCDIR) $(TESTDIR)

EXTRA_DIST = autogen.sh itpp.sln itpp.vcproj Makefile.common TODO

bin_SCRIPTS = @PACKAGE@-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @PACKAGE@.pc

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/@PACKAGE@
