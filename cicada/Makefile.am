AM_CPPFLAGS = -I.. $(BOOST_CPPFLAGS) $(ICU_CPPFLAGS) $(PTHREAD_FLAGS) $(SPARSEHASH_CPPFLAGS) 

LIBUTILS=$(top_builddir)/utils/libutils.la

lib_LTLIBRARIES = libcicada.la

noinst_LTLIBRARIES = \
libcicada-feature.la \
libcicada-semiring.la \
libcicada-eval.la

libcicadaincludedir= $(pkgincludedir)/cicada

semiringdir = $(pkgincludedir)/cicada/semiring
evaldir     = $(pkgincludedir)/cicada/eval
optimizedir = $(pkgincludedir)/cicada/optimize
featuredir  = $(pkgincludedir)/cicada/feature

dist_semiring_HEADERS = \
semiring/envelope.hpp \
semiring/expectation.hpp \
semiring/log.hpp \
semiring/logprob.hpp \
semiring/tropical.hpp \
semiring/traits.hpp

libcicada_semiring_la_SOURCES = \
semiring/envelope.cpp
libcicada_semiring_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_eval_HEADERS = \
eval/bleu.hpp \
eval/score.hpp 

libcicada_eval_la_SOURCES = \
eval/score.cpp
libcicada_eval_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_optimize_HEADERS = \
optimize/line_search.hpp \
optimize/powell.hpp

dist_feature_HEADERS = \
feature/bleu.hpp \
feature/bleu_linear.hpp \
feature/neighbours.hpp \
feature/ngram.hpp \
feature/penalty.hpp \
feature/variational.hpp

libcicada_feature_la_SOURCES = \
feature/bleu.cpp \
feature/bleu_linear.cpp \
feature/neighbours.cpp \
feature/ngram.cpp \
feature/variational.cpp
libcicada_feature_la_CPPFLAGS = $(AM_CPPFLAGS)

libcicadainclude_HEADERS = \
apply.hpp \
apply_cube_prune.hpp \
apply_exact.hpp \
binarize.hpp \
compose.hpp \
compose_cky.hpp \
eval.hpp \
feature.hpp \
feature_function.hpp \
feature_vector.hpp \
ngram.hpp \
grammar.hpp \
grammar_mutable.hpp \
grammar_static.hpp \
hypergraph.hpp \
inside_outside.hpp \
kbest.hpp \
lattice.hpp \
model.hpp \
ngram.hpp \
ngram_index.hpp \
optimize.hpp \
parameter.hpp \
permute.hpp \
prune.hpp \
quantizer.hpp \
rule.hpp \
semiring.hpp \
sentence.hpp \
sort.hpp \
symbol.hpp \
symbol_vector.hpp \
transducer.hpp \
viterbi.hpp \
vocab.hpp \
weight_vector.hpp

libcicada_la_SOURCES = \
feature.cpp \
feature_function.cpp \
grammar_mutable.cpp \
grammar_static.cpp \
hypergraph.cpp \
lattice.cpp \
ngram.cpp \
ngram_index.cpp \
model.cpp \
parameter.cpp \
rule.cpp \
sentence.cpp \
symbol.cpp \
symbol_vector.cpp \
vocab.cpp

libcicada_la_LDFLAGS = \
	$(PTHREAD_FLAGS) \
	$(BOOST_THREAD_LDFLAGS) \
	$(BOOST_FILESYSTEM_LDFLAGS) \
	$(BOOST_IOSTREAMS_LDFLAGS) \
	$(BOOST_SYSTEM_LDFLAGS)

libcicada_la_LIBADD = \
	libcicada-eval.la \
	libcicada-feature.la \
	libcicada-semiring.la \
	\
	$(LIBUTILS) \
	$(BOOST_THREAD_LIBS) \
	$(BOOST_FILESYSTEM_LIBS) \
	$(BOOST_IOSTREAMS_LIBS) \
	$(BOOST_SYSTEM_LIBS) \
	$(ICU_LDFLAGS)

noinst_PROGRAMS = \
grammar_mutable_main \
grammar_static_main \
hypergraph_main \
lattice_main \
parameter_main \
rule_main \
sentence_main \
symbol_main

grammar_mutable_main_SOURCES = grammar_mutable_main.cpp
grammar_mutable_main_LDADD = libcicada.la

grammar_static_main_SOURCES = grammar_static_main.cpp
grammar_static_main_LDADD = libcicada.la

hypergraph_main_SOURCES = hypergraph_main.cpp
hypergraph_main_LDADD = libcicada.la

lattice_main_SOURCES = lattice_main.cpp
lattice_main_LDADD = libcicada.la

parameter_main_SOURCES = parameter_main.cpp
parameter_main_LDADD = libcicada.la

rule_main_SOURCES = rule_main.cpp
rule_main_LDADD = libcicada.la

sentence_main_SOURCES = sentence_main.cpp
sentence_main_LDADD = libcicada.la

symbol_main_SOURCES = symbol_main.cpp
symbol_main_LDADD = libcicada.la
