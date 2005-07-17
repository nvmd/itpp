SUBDIRS = images tutorial local

docdir = $(datadir)/doc/@PACKAGE@-@VERSION@

all-local: html-local

html-local: html/index.html
	if test -d html; then \
		cp $(srcdir)/images/itpp_logga.jpg html; \
	fi

html/index.html: doxygen_html.cfg
	doxygen $<

install-data-local: html-local
	$(mkinstalldirs) $(DESTDIR)$(docdir)
	if test -d html; then \
		$(mkinstalldirs) $(DESTDIR)$(docdir)/html; \
		for f in html/*; do \
			$(INSTALL_DATA) $$f $(DESTDIR)$(docdir)/html; \
		done \
	fi

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)

clean-local:
	rm -rf html *.log *.tag
