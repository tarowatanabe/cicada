AM_CPPFLAGS = -I..  $(ZLIB_CPPFLAGS) $(BZLIB_CPPFLAGS) $(LZMALIB_CPPFLAGS) $(BOOST_CPPFLAGS) $(TCMALLOC_CPPFLAGS) $(JEMALLOC_CPPFLAGS) $(PTHREAD_FLAGS)

LIBKENUTIL=libkenutil.la

noinst_LTLIBRARIES=libkenutil.la

libkenutil_la_SOURCES = \
stream/chain.cc \
stream/io.cc \
stream/line_input.cc \
stream/multi_progress.cc \
double-conversion/bignum-dtoa.cc \
double-conversion/bignum.cc \
double-conversion/cached-powers.cc \
double-conversion/diy-fp.cc \
double-conversion/double-conversion.cc \
double-conversion/fast-dtoa.cc \
double-conversion/fixed-dtoa.cc \
double-conversion/strtod.cc \
bit_packing.cc \
ersatz_progress.cc \
exception.cc \
file.cc \
file_piece.cc \
mmap.cc \
murmur_hash.cc \
pool.cc \
read_compressed.cc \
scoped.cc \
string_piece.cc \
usage.cc

libkenutil_la_LDFLAGS = -version-info $(CICADA_LTVERSION)

libkenutil_la_LIBADD = \
	$(ZLIB_LDFLAGS) \
	$(BZLIB_LDFLAGS) \
	$(LZMALIB_LDFLAGS) \
	$(BOOST_THREAD_LDFLAGS) \
	$(JEMALLOC_LDFLAGS) \
	$(TCMALLOC_LDFLAGS)

noinst_HEADERS = \
stream/block.hh \
stream/chain.hh \
stream/config.hh \
stream/io.hh \
stream/line_input.hh \
stream/multi_progress.hh \
stream/sort.hh \
stream/stream.hh \
stream/timer.hh \
double-conversion/bignum-dtoa.h \
double-conversion/bignum.h \
double-conversion/cached-powers.h \
double-conversion/diy-fp.h \
double-conversion/double-conversion.h \
double-conversion/fast-dtoa.h \
double-conversion/fixed-dtoa.h \
double-conversion/ieee.h \
double-conversion/strtod.h \
double-conversion/utils.h \
bit_packing.hh \
ersatz_progress.hh \
exception.hh \
fake_ofstream.hh \
file.hh \
file_piece.hh \
getopt.hh \
have.hh \
joint_sort.hh \
mmap.hh \
multi_intersection.hh \
murmur_hash.hh \
pcqueue.hh \
pool.hh \
probing_hash_table.hh \
proxy_iterator.hh \
read_compressed.hh \
scoped.hh \
sized_iterator.hh \
sorted_uniform.hh \
string_piece.hh \
string_piece_hash.hh \
thread_pool.hh \
tokenize_piece.hh \
usage.hh

EXTRA_DIST = \
Jamfile
