AM_CPPFLAGS = -I.. $(KENLM_CPPFLAGS) $(BOOST_CPPFLAGS) $(TCMALLOC_CPPFLAGS) $(JEMALLOC_CPPFLAGS) $(PTHREAD_FLAGS)

LIBKENLM=libkenlm.la
LIBKENUTIL=../util/libkenutil.la

noinst_LTLIBRARIES=libkenlm.la

libkenlm_la_SOURCES = \
bhiksha.cc \
binary_format.cc \
config.cc \
lm_exception.cc \
model.cc \
quantize.cc \
read_arpa.cc \
search_hashed.cc \
search_trie.cc \
sizes.cc \
trie.cc \
trie_sort.cc \
value_build.cc \
virtual_interface.cc \
vocab.cc

libkenlm_la_LDFLAGS = -version-info $(CICADA_LTVERSION)

libkenlm_la_LIBADD = \
	$(JEMALLOC_LDFLAGS) \
	$(TCMALLOC_LDFLAGS) \
	$(LIBKENUTIL)

noinst_HEADERS = \
bhiksha.hh \
binary_format.hh \
blank.hh \
config.hh \
enumerate_vocab.hh \
facade.hh \
left.hh \
lm_exception.hh \
max_order.hh \
model.hh \
model_type.hh \
ngram_query.hh \
partial.hh \
quantize.hh \
read_arpa.hh \
return.hh \
search_hashed.hh \
search_trie.hh \
sizes.hh \
state.hh \
trie.hh \
trie_sort.hh \
value.hh \
value_build.hh \
virtual_interface.hh \
vocab.hh \
weights.hh \
word_index.hh

