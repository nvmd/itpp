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

dist-html:
	$(MAKE) -C doc dist-html

upload: upload-src upload-html

upload-src: dist
	@UPLOAD_COMMAND@ @PACKAGE_TARNAME@-@VERSION@.tar.gz \
		@PACKAGE_TARNAME@-@VERSION@.tar.bz2

upload-html: dist-html
	@UPLOAD_COMMAND@ @PACKAGE_TARNAME@-html-doc-@VERSION@.tar.gz \
		@PACKAGE_TARNAME@-html-doc-@VERSION@.tar.bz2
