AM_CPPFLAGS = -I.. $(BOOST_CPPFLAGS) $(SNAPPY_CPPFLAGS)

LIBCODEC=libcodec.la

lib_LTLIBRARIES = libcodec.la

libcodec_la_SOURCES = \
fastlz.c \
lz4.c \
lz4hc.c \
lz4_codec.cpp \
quicklz.c \
quicklz_codec.cpp

libcodec_la_LDFLAGS = -version-info $(CICADA_LTVERSION)

libcodec_la_LIBADD = \
	$(BOOST_IOSTREAMS_LDFLAGS) $(BOOST_IOSTREAMS_LIBS) \
	$(SNAPPY_LDFLAGS)

libcodecincludedir = $(pkgincludedir)/codec
libcodecinclude_HEADERS = \
lz4.hpp

if WITH_SNAPPY
libcodec_la_SOURCES += snappy_codec.cpp
libcodecinclude_HEADERS += snappy.hpp
endif WITH_SNAPPY

noinst_HEADERS = \
codec_impl.hpp \
fastlz.h \
lz4.h \
lz4hc.h \
quicklz.h

noinst_PROGRAMS = \
lz4_compress_main \
lz4_decompress_main \
quicklz_compress_main \
quicklz_decompress_main \
snappy_compress_main \
snappy_decompress_main

lz4_compress_main_SOURCES = lz4_compress_main.cpp
lz4_compress_main_LDADD = $(LIBCODEC)

lz4_decompress_main_SOURCES = lz4_decompress_main.cpp
lz4_decompress_main_LDADD = $(LIBCODEC)

quicklz_compress_main_SOURCES = quicklz_compress_main.cpp
quicklz_compress_main_LDADD = $(LIBCODEC)

quicklz_decompress_main_SOURCES = quicklz_decompress_main.cpp
quicklz_decompress_main_LDADD = $(LIBCODEC)

snappy_compress_main_SOURCES = snappy_compress_main.cpp
snappy_compress_main_LDADD = $(LIBCODEC)

snappy_decompress_main_SOURCES = snappy_decompress_main.cpp
snappy_decompress_main_LDADD = $(LIBCODEC)

noinst_DATA = .gitignore

.gitignore: Makefile.am
	cat $(top_builddir)/.gitignore > .gitignore
	for bin in $(PROGRAMS); do \
	  echo $$bin >> .gitignore; \
	done; \
	for bin in $(SCRIPTS); do \
	  echo $$bin >> .gitignore; \
	done; \
	for bin in $(IGNORES); do \
	  echo $$bin >> .gitignore; \
	done
