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

EXTRA_DIST = autogen.sh ChangeLog-2005 ChangeLog-2006 configure.ac.in \
	itpp-config.1.in itpp.spec itpp.spec.in Makefile.common NEWS-3.10 \
	NEWS-3.9 NEWS-3.99 TODO VERSION

man_MANS = itpp-config.1
dist_man_MANS = itpp-config.1

bin_SCRIPTS = itpp-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = itpp.pc

clean-local:
	rm -f config_summary.log

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/itpp

dist-html:
	$(MAKE) -C doc dist-html

snapshot_suffix = _pre`date +"%Y%m%d"`

snapshot:
	$(MAKE) dist distdir=@PACKAGE@-@VERSION@$(snapshot_suffix)

snapshot-html:
	$(MAKE) -C doc dist-html \
		distdir_html=@PACKAGE@-html-doc-@VERSION@$(snapshot_suffix)
