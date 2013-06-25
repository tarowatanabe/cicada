AM_CPPFLAGS = -I.. $(KENLM_CPPFLAGS) $(BOOST_CPPFLAGS) $(TCMALLOC_CPPFLAGS) $(JEMALLOC_CPPFLAGS) $(PTHREAD_FLAGS)

LIBKENLM=../lm/libkenlm.la
LIBKENUTIL=../util/libkenutil.la

noinst_LTLIBRARIES=libkensearch.la

libkensearch_la_SOURCES = \
edge_generator.cc \
nbest.cc \
rule.cc \
vertex.cc

noinst_HEADERS = \
applied.hh \
config.hh \
context.hh \
edge.hh \
edge_generator.hh \
header.hh \
nbest.hh \
rule.hh \
types.hh \
vertex.hh \
vertex_generator.hh

libkensearch_la_LDFLAGS = -version-info $(CICADA_LTVERSION)

libkensearch_la_LIBADD = \
	$(LIBKENLM)
