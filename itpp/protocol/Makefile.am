include $(top_srcdir)/Makefile.common

h_sources = \
	events.h \
	front_drop_queue.h \
	packet_channel.h \
	packet_generator.h \
	packet.h \
	selective_repeat.h \
	signals_slots.h \
	tcp_client_server.h \
	tcp.h

cpp_sources = \
	events.cpp \
	front_drop_queue.cpp \
	packet_channel.cpp \
	packet_generator.cpp \
	selective_repeat.cpp 
	tcp.cpp

noinst_LTLIBRARIES = libprotocol.la
libprotocol_la_SOURCES = $(h_sources) $(cpp_sources) 

pkgincludedir = $(includedir)/@PACKAGE@/protocol
pkginclude_HEADERS = $(h_sources)
