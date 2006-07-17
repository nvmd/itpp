include $(top_srcdir)/Makefile.common

ACLOCAL_AMFLAGS = -I config

AUTOMAKE_OPTIONS = dist-bzip2

if HTML_DOC
  HTML_DOC_DIR = doc
endif

if RUN_TESTS
  TEST_DIR = tests
endif

SUBDIRS = itpp extras win32 $(HTML_DOC_DIR) $(TEST_DIR)

EXTRA_DIST = autogen.sh itpp.spec Makefile.common TODO

bin_SCRIPTS = @PACKAGE@-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @PACKAGE@.pc

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/@PACKAGE@

disthtml:
	$(MAKE) -C doc disthtml
