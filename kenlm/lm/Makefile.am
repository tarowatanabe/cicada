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

bin_PROGRAMS = \
cicada_kenlm_build_binary \
cicada_kenlm_fragment \
cicada_kenlm_query \
cicada_kenlm_lmplz

cicada_kenlm_build_binary_SOURCES = build_binary_main.cc
cicada_kenlm_build_binary_LDADD   = $(LIBKENLM)

cicada_kenlm_fragment_SOURCES = fragment_main.cc
cicada_kenlm_fragment_LDADD   = $(LIBKENLM)

cicada_kenlm_query_SOURCES = query_main.cc
cicada_kenlm_query_LDADD   = $(LIBKENLM)

cicada_kenlm_lmplz_SOURCES = \
builder/adjust_counts.hh \
builder/corpus_count.hh \
builder/discount.hh \
builder/header_info.hh \
builder/initial_probabilities.hh \
builder/interpolate.hh \
builder/joint_order.hh \
builder/multi_stream.hh \
builder/ngram.hh \
builder/ngram_stream.hh \
builder/pipeline.hh \
builder/print.hh \
builder/sort.hh \
builder/adjust_counts.cc \
builder/corpus_count.cc \
builder/initial_probabilities.cc \
builder/interpolate.cc \
builder/pipeline.cc \
builder/print.cc \
builder/lmplz_main.cc

cicada_kenlm_lmplz_CPPFLAGS = $(AM_CPPFLAGS)

cicada_kenlm_lmplz_LDADD = \
$(LIBKENLM) \
$(BOOST_THREAD_LDFLAGS) \
$(BOOST_THREAD_LIBS) \
$(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
$(BOOST_PROGRAM_OPTIONS_LIBS)

EXTRA_DIST = \
Jamfile

noinst_DATA = .gitignore

.gitignore: Makefile.am
	cat $(top_builddir)/.gitignore > .gitignore
	for bin in $(PROGRAMS); do \
	  echo $$bin >> .gitignore; \
	done; \
	for bin in $(SCRIPTS); do \
	  echo $$bin >> .gitignore; \
	done
