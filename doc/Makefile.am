include $(top_srcdir)/doc/local/sources.mk
include $(top_srcdir)/doc/tutorial/sources.mk
include $(top_srcdir)/doc/tutorial/src/sources.mk
include $(top_srcdir)/itpp/sources.mk
include $(top_srcdir)/itpp/base/sources.mk
include $(top_srcdir)/itpp/base/bessel/sources.mk
include $(top_srcdir)/itpp/comm/sources.mk
include $(top_srcdir)/itpp/fixedpoint/sources.mk
include $(top_srcdir)/itpp/protocol/sources.mk
include $(top_srcdir)/itpp/srccode/sources.mk

SUBDIRS = images tutorial local

docdir = @docdir@

all-local: html-local

html-local: html/index.html

html/index.html: doxygen_html.cfg \
	$(doc_local_sources) $(html_local_sources) \
	$(doc_tutorial_sources) $(cpp_tutorial_sources) \
	$(h_sources) $(cpp_sources) \
	$(h_base_sources) $(cpp_base_sources) \
	$(h_base_bessel_sources) $(cpp_base_bessel_sources) \
	$(h_comm_sources) $(cpp_comm_sources) \
	$(h_fixedpoint_sources)	$(cpp_fixedpoint_sources) \
	$(h_protocol_sources) $(cpp_protocol_sources) \
	$(h_srccode_sources) $(cpp_srccode_sources)
	doxygen $<; \
	if test -d html; then \
		cp $(srcdir)/images/itpp_logga.jpg html; \
	fi

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
