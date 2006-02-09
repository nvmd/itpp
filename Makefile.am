include $(top_srcdir)/Makefile.common

ACLOCAL_AMFLAGS = -I config

if HTML_DOC
HTMLDOCDIR = doc
endif

if HAVE_DIFF
TESTDIR = tests
endif

SUBDIRS = itpp extras $(HTMLDOCDIR) $(TESTDIR)

EXTRA_DIST = autogen.sh itpp.spec Makefile.common TODO \
	itpp.sln itpp.vcproj itpp_msvc2003.sln itpp_msvc2003.vcproj

bin_SCRIPTS = @PACKAGE@-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @PACKAGE@.pc

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/@PACKAGE@

disthtml:
	$(MAKE) -C doc disthtml
