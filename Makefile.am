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

EXTRA_DIST = autogen.sh ChangeLog-2005 ChangeLog-2006 itpp.spec \
	Makefile.common NEWS-3.9 NEWS-3.10 TODO

man_MANS = @PACKAGE@-config.1
dist_man_MANS = @PACKAGE@-config.1

bin_SCRIPTS = @PACKAGE@-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @PACKAGE@.pc

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/@PACKAGE@

dist-html:
	$(MAKE) -C doc dist-html


snapshot: snapshot-src snapshot-html

snapshot-src:
	$(MAKE) dist distdir=@PACKAGE_TARNAME@-`date +"%Y%m%d"`

snapshot-html:
	$(MAKE) -C doc dist-html html_tarname=@PACKAGE_TARNAME@-html-doc-`date +"%Y%m%d"`


upload: upload-src upload-html

upload-src: dist
	@UPLOAD_COMMAND@ @PACKAGE_TARNAME@-@VERSION@.tar.gz \
		@PACKAGE_TARNAME@-@VERSION@.tar.bz2

upload-html: dist-html
	@UPLOAD_COMMAND@ @PACKAGE_TARNAME@-html-doc-@VERSION@.tar.gz \
		@PACKAGE_TARNAME@-html-doc-@VERSION@.tar.bz2
