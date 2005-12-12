SUBDIRS = images tutorial local

docdir = @docdir@

all-local: html-local

html-local: html/index.html

html/index.html: doxygen_html.cfg \
	$(doc_local_sources) \
	$(doc_tutorial_sources) $(cpp_tutorial_sources) \
	$(h_sources) \
	$(h_base_sources) $(cpp_base_sources) \
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

doc_local_sources = \
	$(srcdir)/local/codingrules.doc \
	$(srcdir)/local/copyright.doc \
	$(srcdir)/local/documentation.doc \
	$(srcdir)/local/help_wanted.doc \
	$(srcdir)/local/index.doc \
	$(srcdir)/local/installation.doc \
	$(srcdir)/local/linking.doc \
	$(srcdir)/local/test.doc \
	$(srcdir)/local/users_guide.doc \
	$(srcdir)/local/verification.doc

doc_tutorial_sources = \
	$(srcdir)/tutorial/bpsk.doc \
	$(srcdir)/tutorial/convcode.doc \
	$(srcdir)/tutorial/interleaver.doc \
	$(srcdir)/tutorial/itfile.doc \
	$(srcdir)/tutorial/qpsk_simulation.doc \
	$(srcdir)/tutorial/rayleigh.doc \
	$(srcdir)/tutorial/reedsolomon.doc \
	$(srcdir)/tutorial/spread.doc \
	$(srcdir)/tutorial/timer.doc \
	$(srcdir)/tutorial/tutorial.doc \
	$(srcdir)/tutorial/vector_and_matrix.doc

h_base_sources = \
	$(top_srcdir)/itpp/base/array.h \
	$(top_srcdir)/itpp/base/bessel/bessel_internal.h \
	$(top_srcdir)/itpp/base/bessel.h \
	$(top_srcdir)/itpp/base/binary.h \
	$(top_srcdir)/itpp/base/binfile.h \
	$(top_srcdir)/itpp/base/cblas.h \
	$(top_srcdir)/itpp/base/cholesky.h \
	$(top_srcdir)/itpp/base/circular_buffer.h \
	$(top_srcdir)/itpp/base/converters.h \
	$(top_srcdir)/itpp/base/copy_vector.h \
	$(top_srcdir)/itpp/base/det.h \
	$(top_srcdir)/itpp/base/eigen.h \
	$(top_srcdir)/itpp/base/elmatfunc.h \
	$(top_srcdir)/itpp/base/factory.h \
	$(top_srcdir)/itpp/base/fastica.h \
	$(top_srcdir)/itpp/base/fastmath.h \
	$(top_srcdir)/itpp/base/fftw3.h \
	$(top_srcdir)/itpp/base/filter.h \
	$(top_srcdir)/itpp/base/filter_design.h \
	$(top_srcdir)/itpp/base/freq_filt.h \
	$(top_srcdir)/itpp/base/gf2mat.h \
	$(top_srcdir)/itpp/base/help_functions.h \
	$(top_srcdir)/itpp/base/integration.h \
	$(top_srcdir)/itpp/base/inv.h \
	$(top_srcdir)/itpp/base/itassert.h \
	$(top_srcdir)/itpp/base/itfile.h \
	$(top_srcdir)/itpp/base/itpp_version.h \
	$(top_srcdir)/itpp/base/lapack.h \
	$(top_srcdir)/itpp/base/ls_solve.h \
	$(top_srcdir)/itpp/base/lu.h \
	$(top_srcdir)/itpp/base/machdep.h \
	$(top_srcdir)/itpp/base/matfunc.h \
	$(top_srcdir)/itpp/base/mat.h \
	$(top_srcdir)/itpp/base/operators.h \
	$(top_srcdir)/itpp/base/parser.h \
	$(top_srcdir)/itpp/base/poly.h \
	$(top_srcdir)/itpp/base/qr.h \
	$(top_srcdir)/itpp/base/random.h \
	$(top_srcdir)/itpp/base/scalfunc.h \
	$(top_srcdir)/itpp/base/sigfun.h \
	$(top_srcdir)/itpp/base/smat.h \
	$(top_srcdir)/itpp/base/sort.h \
	$(top_srcdir)/itpp/base/source.h \
	$(top_srcdir)/itpp/base/specmat.h \
	$(top_srcdir)/itpp/base/stack.h \
	$(top_srcdir)/itpp/base/stat.h \
	$(top_srcdir)/itpp/base/svd.h \
	$(top_srcdir)/itpp/base/svec.h \
	$(top_srcdir)/itpp/base/timing.h \
	$(top_srcdir)/itpp/base/transforms.h \
	$(top_srcdir)/itpp/base/vec.h

h_comm_sources = \
	$(top_srcdir)/itpp/comm/bch.h \
	$(top_srcdir)/itpp/comm/channel_code.h \
	$(top_srcdir)/itpp/comm/channel.h \
	$(top_srcdir)/itpp/comm/commfunc.h \
	$(top_srcdir)/itpp/comm/convcode.h \
	$(top_srcdir)/itpp/comm/crc.h \
	$(top_srcdir)/itpp/comm/egolay.h \
	$(top_srcdir)/itpp/comm/error_counters.h \
	$(top_srcdir)/itpp/comm/galois.h \
	$(top_srcdir)/itpp/comm/hammcode.h \
	$(top_srcdir)/itpp/comm/interleave.h \
	$(top_srcdir)/itpp/comm/modulator.h \
	$(top_srcdir)/itpp/comm/ofdm.h \
	$(top_srcdir)/itpp/comm/pulse_shape.h \
	$(top_srcdir)/itpp/comm/punct_convcode.h \
	$(top_srcdir)/itpp/comm/rec_syst_conv_code.h \
	$(top_srcdir)/itpp/comm/reedsolomon.h \
	$(top_srcdir)/itpp/comm/sequence.h \
	$(top_srcdir)/itpp/comm/spread.h \
	$(top_srcdir)/itpp/comm/turbo.h

h_fixedpoint_sources = \
	$(top_srcdir)/itpp/fixedpoint/cfixed.h \
	$(top_srcdir)/itpp/fixedpoint/cfix.h \
	$(top_srcdir)/itpp/fixedpoint/fix_base.h \
	$(top_srcdir)/itpp/fixedpoint/fixed.h \
	$(top_srcdir)/itpp/fixedpoint/fix_factory.h \
	$(top_srcdir)/itpp/fixedpoint/fix_functions.h \
	$(top_srcdir)/itpp/fixedpoint/fix.h \
	$(top_srcdir)/itpp/fixedpoint/fix_operators.h

h_protocol_sources = \
	$(top_srcdir)/itpp/protocol/events.h \
	$(top_srcdir)/itpp/protocol/front_drop_queue.h \
	$(top_srcdir)/itpp/protocol/packet_channel.h \
	$(top_srcdir)/itpp/protocol/packet_generator.h \
	$(top_srcdir)/itpp/protocol/packet.h \
	$(top_srcdir)/itpp/protocol/selective_repeat.h \
	$(top_srcdir)/itpp/protocol/signals_slots.h \
	$(top_srcdir)/itpp/protocol/tcp_client_server.h \
	$(top_srcdir)/itpp/protocol/tcp.h

h_srccode_sources = \
	$(top_srcdir)/itpp/srccode/audiofile.h \
	$(top_srcdir)/itpp/srccode/gmm.h \
	$(top_srcdir)/itpp/srccode/lpcfunc.h \
	$(top_srcdir)/itpp/srccode/pnm.h \
	$(top_srcdir)/itpp/srccode/vq.h \
	$(top_srcdir)/itpp/srccode/vqtrain.h

h_sources = \
	$(top_srcdir)/itpp/config.h \
	$(top_srcdir)/itpp/itbase.h \
	$(top_srcdir)/itpp/itcomm.h \
	$(top_srcdir)/itpp/itconfig.h \
	$(top_srcdir)/itpp/itfixedpoint.h \
	$(top_srcdir)/itpp/itmex.h \
	$(top_srcdir)/itpp/itprotocol.h \
	$(top_srcdir)/itpp/itsrccode.h

cpp_base_sources = \
	$(top_srcdir)/itpp/base/bessel/airy.cpp \
	$(top_srcdir)/itpp/base/bessel/chbevl.cpp \
	$(top_srcdir)/itpp/base/bessel.cpp \
	$(top_srcdir)/itpp/base/bessel/hyperg.cpp \
	$(top_srcdir)/itpp/base/bessel/i0.cpp \
	$(top_srcdir)/itpp/base/bessel/i1.cpp \
	$(top_srcdir)/itpp/base/bessel/iv.cpp \
	$(top_srcdir)/itpp/base/bessel/jv.cpp \
	$(top_srcdir)/itpp/base/bessel/k0.cpp \
	$(top_srcdir)/itpp/base/bessel/k1.cpp \
	$(top_srcdir)/itpp/base/bessel/kn.cpp \
	$(top_srcdir)/itpp/base/bessel/polevl.cpp \
	$(top_srcdir)/itpp/base/bessel/struve.cpp \
	$(top_srcdir)/itpp/base/binary.cpp \
	$(top_srcdir)/itpp/base/binfile.cpp \
	$(top_srcdir)/itpp/base/cholesky.cpp \
	$(top_srcdir)/itpp/base/converters.cpp \
	$(top_srcdir)/itpp/base/det.cpp \
	$(top_srcdir)/itpp/base/eigen.cpp \
	$(top_srcdir)/itpp/base/elmatfunc.cpp \
	$(top_srcdir)/itpp/base/fastica.cpp \
	$(top_srcdir)/itpp/base/fastmath.cpp \
	$(top_srcdir)/itpp/base/filter.cpp \
	$(top_srcdir)/itpp/base/filter_design.cpp \
	$(top_srcdir)/itpp/base/freq_filt.cpp \
	$(top_srcdir)/itpp/base/gf2mat.cpp \
	$(top_srcdir)/itpp/base/integration.cpp \
	$(top_srcdir)/itpp/base/inv.cpp \
	$(top_srcdir)/itpp/base/itassert.cpp \
	$(top_srcdir)/itpp/base/itfile.cpp \
	$(top_srcdir)/itpp/base/itpp_version.cpp \
	$(top_srcdir)/itpp/base/ls_solve.cpp \
	$(top_srcdir)/itpp/base/lu.cpp \
	$(top_srcdir)/itpp/base/mat.cpp \
	$(top_srcdir)/itpp/base/matfunc.cpp \
	$(top_srcdir)/itpp/base/operators.cpp \
	$(top_srcdir)/itpp/base/parser.cpp \
	$(top_srcdir)/itpp/base/poly.cpp \
	$(top_srcdir)/itpp/base/qr.cpp \
	$(top_srcdir)/itpp/base/random.cpp \
	$(top_srcdir)/itpp/base/scalfunc.cpp \
	$(top_srcdir)/itpp/base/sigfun.cpp \
	$(top_srcdir)/itpp/base/smat.cpp \
	$(top_srcdir)/itpp/base/source.cpp \
	$(top_srcdir)/itpp/base/specmat.cpp \
	$(top_srcdir)/itpp/base/stat.cpp \
	$(top_srcdir)/itpp/base/svd.cpp \
	$(top_srcdir)/itpp/base/svec.cpp \
	$(top_srcdir)/itpp/base/timing.cpp \
	$(top_srcdir)/itpp/base/transforms.cpp \
	$(top_srcdir)/itpp/base/vec.cpp

cpp_comm_sources = \
	$(top_srcdir)/itpp/comm/bch.cpp \
	$(top_srcdir)/itpp/comm/channel.cpp \
	$(top_srcdir)/itpp/comm/commfunc.cpp \
	$(top_srcdir)/itpp/comm/convcode.cpp \
	$(top_srcdir)/itpp/comm/crc.cpp \
	$(top_srcdir)/itpp/comm/egolay.cpp \
	$(top_srcdir)/itpp/comm/error_counters.cpp \
	$(top_srcdir)/itpp/comm/galois.cpp \
	$(top_srcdir)/itpp/comm/hammcode.cpp \
	$(top_srcdir)/itpp/comm/interleave.cpp \
	$(top_srcdir)/itpp/comm/modulator.cpp \
	$(top_srcdir)/itpp/comm/ofdm.cpp \
	$(top_srcdir)/itpp/comm/pulse_shape.cpp \
	$(top_srcdir)/itpp/comm/punct_convcode.cpp \
	$(top_srcdir)/itpp/comm/rec_syst_conv_code.cpp \
	$(top_srcdir)/itpp/comm/reedsolomon.cpp \
	$(top_srcdir)/itpp/comm/sequence.cpp \
	$(top_srcdir)/itpp/comm/spread.cpp \
	$(top_srcdir)/itpp/comm/turbo.cpp

cpp_fixedpoint_sources = \
	$(top_srcdir)/itpp/fixedpoint/cfix.cpp \
	$(top_srcdir)/itpp/fixedpoint/cfixed.cpp \
	$(top_srcdir)/itpp/fixedpoint/fix_base.cpp \
	$(top_srcdir)/itpp/fixedpoint/fix.cpp \
	$(top_srcdir)/itpp/fixedpoint/fixed.cpp \
	$(top_srcdir)/itpp/fixedpoint/fix_factory.cpp \
	$(top_srcdir)/itpp/fixedpoint/fix_functions.cpp \
	$(top_srcdir)/itpp/fixedpoint/fix_operators.cpp

cpp_protocol_sources = \
	$(top_srcdir)/itpp/protocol/events.cpp \
	$(top_srcdir)/itpp/protocol/front_drop_queue.cpp \
	$(top_srcdir)/itpp/protocol/packet_channel.cpp \
	$(top_srcdir)/itpp/protocol/packet_generator.cpp \
	$(top_srcdir)/itpp/protocol/selective_repeat.cpp \
	$(top_srcdir)/itpp/protocol/tcp.cpp

cpp_srccode_sources = \
	$(top_srcdir)/itpp/srccode/audiofile.cpp \
	$(top_srcdir)/itpp/srccode/gmm.cpp \
	$(top_srcdir)/itpp/srccode/lpcfunc.cpp \
	$(top_srcdir)/itpp/srccode/pnm.cpp \
	$(top_srcdir)/itpp/srccode/vq.cpp \
	$(top_srcdir)/itpp/srccode/vqtrain.cpp

cpp_tutorial_sources = \
	$(srcdir)/tutorial/src/bpsk.cpp \
	$(srcdir)/tutorial/src/convcode.cpp \
	$(srcdir)/tutorial/src/interleaver.cpp \
	$(srcdir)/tutorial/src/qpsk_simulation.cpp \
	$(srcdir)/tutorial/src/rayleigh.cpp \
	$(srcdir)/tutorial/src/read_it_file.cpp \
	$(srcdir)/tutorial/src/reedsolomon.cpp \
	$(srcdir)/tutorial/src/spread.cpp \
	$(srcdir)/tutorial/src/timer.cpp \
	$(srcdir)/tutorial/src/vector_and_matrix.cpp \
	$(srcdir)/tutorial/src/write_it_file.cpp
