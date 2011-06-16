AM_CPPFLAGS = -I.. $(BOOST_CPPFLAGS) $(ICU_CPPFLAGS) $(PTHREAD_FLAGS) $(SPARSEHASH_CPPFLAGS)

LIBUTILS=$(top_builddir)/utils/libutils.la

lib_LTLIBRARIES = libcicada.la

noinst_LTLIBRARIES = \
libcicada-feature.la \
libcicada-matcher.la \
libcicada-semiring.la \
libcicada-eval.la \
libcicada-format.la \
libcicada-head.la \
libcicada-optimize.la \
libcicada-operation.la \
libcicada-signature.la  \
libcicada-stemmer.la 

libcicadaincludedir= $(pkgincludedir)/cicada

semiringdir  = $(pkgincludedir)/cicada/semiring
evaldir      = $(pkgincludedir)/cicada/eval
formatdir    = $(pkgincludedir)/cicada/format
headdir      = $(pkgincludedir)/cicada/head
optimizedir  = $(pkgincludedir)/cicada/optimize
featuredir   = $(pkgincludedir)/cicada/feature
signaturedir = $(pkgincludedir)/cicada/signature
stemmerdir   = $(pkgincludedir)/cicada/stemmer
operationdir = $(pkgincludedir)/cicada/operation
tokenizerdir = $(pkgincludedir)/cicada/tokenizer
matcherdir   = $(pkgincludedir)/cicada/matcher

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
eval/combined.hpp \
eval/f.hpp \
eval/bleu.hpp \
eval/parseval.hpp \
eval/per.hpp \
eval/ribes.hpp \
eval/score.hpp \
eval/sb.hpp \
eval/sk.hpp \
eval/ter.hpp \
eval/wer.hpp \
eval/wlcs.hpp

libcicada_eval_la_SOURCES = \
eval/bleu.cpp \
eval/combined.cpp \
eval/decode.hpp \
eval/encode.hpp \
eval/f.cpp \
eval/parseval.cpp \
eval/per.cpp \
eval/ribes.cpp \
eval/score.cpp \
eval/sb.cpp \
eval/sk.cpp \
eval/ter.cpp \
eval/wer.cpp \
eval/wlcs.cpp
libcicada_eval_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_format_HEADERS = \
format/date.hpp \
format/number.hpp

libcicada_format_la_SOURCES = \
format/date.cpp \
format/number.cpp

libcicada_format_la_LIBADD = $(ICU_LDFLAGS)
libcicada_format_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_head_HEADERS = \
head/chinese.hpp \
head/collins.hpp \
head/collins_modified.hpp 

libcicada_head_la_SOURCES = \
head/chinese.cpp \
head/collins.cpp \
head/collins_modified.cpp

libcicada_head_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_optimize_HEADERS = \
optimize/line_search.hpp \
optimize/powell.hpp

libcicada_optimize_la_SOURCES = \
optimize/line_search.cpp
libcicada_optimize_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_feature_HEADERS = \
feature/alignment.hpp \
feature/antecedent.hpp \
feature/bleu.hpp \
feature/bleu_expected.hpp \
feature/bleu_linear.hpp \
feature/bleu_multi.hpp \
feature/distortion.hpp \
feature/global_lexicon.hpp \
feature/lexicalized_reordering.hpp \
feature/neighbours.hpp \
feature/ngram.hpp \
feature/ngram_tree.hpp \
feature/parent.hpp \
feature/penalty.hpp \
feature/permute.hpp \
feature/rule_shape.hpp \
feature/sgml_tag.hpp \
feature/span.hpp \
feature/variational.hpp \
feature/vocabulary.hpp

libcicada_feature_la_SOURCES = \
feature/alignment.cpp \
feature/antecedent.cpp \
feature/bleu.cpp \
feature/bleu_expected.cpp \
feature/bleu_linear.cpp \
feature/bleu_multi.cpp \
feature/distortion.cpp \
feature/global_lexicon.cpp \
feature/lexicalized_reordering.cpp \
feature/neighbours.cpp \
feature/ngram.cpp \
feature/ngram_tree.cpp \
feature/parent.cpp \
feature/penalty.cpp \
feature/permute.cpp \
feature/rule_shape.cpp \
feature/sgml_tag.cpp \
feature/span.cpp \
feature/variational.cpp \
feature/vocabulary.cpp

libcicada_feature_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_signature_HEADERS = \
signature/chinese.hpp  \
signature/english.hpp 

libcicada_signature_la_SOURCES = \
signature/chinese.cpp \
signature/english.cpp

libcicada_signature_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_stemmer_HEADERS = \
stemmer/arabic.hpp \
stemmer/prefix.hpp \
stemmer/suffix.hpp \
stemmer/latin.hpp \
stemmer/lower.hpp \
stemmer/digit.hpp \
stemmer/nfkc.hpp \
stemmer/snowball.hpp

libcicada_stemmer_la_SOURCES = \
stemmer/arabic.cpp \
stemmer/prefix.cpp \
stemmer/suffix.cpp \
stemmer/latin.cpp \
stemmer/lower.cpp \
stemmer/digit.cpp \
stemmer/nfkc.cpp \
stemmer/snowball.cpp

libcicada_stemmer_la_CPPFLAGS = $(AM_CPPFLAGS)
libcicada_stemmer_la_LIBADD = $(top_builddir)/libstemmer_c/libstemmer.la

dist_tokenizer_HEADERS = \
tokenizer/characters.hpp \
tokenizer/lower.hpp \
tokenizer/nist.hpp \
tokenizer/nonascii.hpp \
tokenizer/penntreebank.hpp \
tokenizer/stemmer.hpp \
tokenizer/tokenize.hpp

dist_matcher_HEADERS = \
matcher/lower.hpp \
matcher/stemmer.hpp \
matcher/wordnet.hpp 

libcicada_matcher_la_SOURCES = \
matcher/wordnet.cpp

libcicada_matcher_la_CPPFLAGS = $(AM_CPPFLAGS)
libcicada_matcher_la_LIBADD = $(top_builddir)/wn/libwn.la

dist_operation_HEADERS = \
operation/apply.hpp \
operation/binarize.hpp \
operation/clear.hpp \
operation/compose.hpp \
operation/debinarize.hpp \
operation/expected_ngram.hpp \
operation/functional.hpp \
operation/generate.hpp \
operation/intersect.hpp \
operation/normalize.hpp \
operation/output.hpp \
operation/parse.hpp \
operation/permute.hpp \
operation/prune.hpp \
operation/push_bos_eos.hpp \
operation/remove_annotation.hpp \
operation/remove_epsilon.hpp \
operation/remove_features.hpp \
operation/sort_tail.hpp \
operation/span_forest.hpp \
operation/traversal.hpp \
operation/viterbi.hpp

libcicada_operation_la_SOURCES = \
operation/apply.cpp \
operation/binarize.cpp \
operation/clear.cpp \
operation/compose.cpp \
operation/debinarize.cpp \
operation/expected_ngram.cpp \
operation/intersect.cpp \
operation/generate.cpp \
operation/normalize.cpp \
operation/output.cpp \
operation/parse.cpp \
operation/permute.cpp \
operation/prune.cpp \
operation/push_bos_eos.cpp \
operation/remove_annotation.cpp \
operation/remove_epsilon.cpp \
operation/remove_feature.cpp \
operation/sort_tail.cpp \
operation/span_forest.cpp \
operation/viterbi.cpp
libcicada_operation_la_CPPFLAGS = $(AM_CPPFLAGS)


libcicadainclude_HEADERS = \
alignment.hpp \
apply.hpp \
apply_cube_grow.hpp \
apply_cube_prune.hpp \
apply_exact.hpp \
apply_incremental.hpp \
apply_state_less.hpp \
attribute.hpp \
attribute_vector.hpp \
binarize.hpp \
binarize_all.hpp \
binarize_base.hpp \
binarize_cyk.hpp \
binarize_left.hpp \
binarize_right.hpp \
binarize_terminal.hpp \
cluster.hpp \
cluster_stemmer.hpp \
compose.hpp \
compose_alignment.hpp \
compose_cky.hpp \
compose_earley.hpp \
compose_grammar.hpp \
compose_phrase.hpp \
compose_tree.hpp \
compose_tree_cky.hpp \
debinarize.hpp \
dot_product.hpp \
edit_distance_forest.hpp \
eval.hpp \
expected_ngram.hpp \
feature.hpp \
feature_function.hpp \
feature_vector.hpp \
format.hpp \
generate.hpp \
generate_earley.hpp \
global_lexicon.hpp \
grammar.hpp \
grammar_format.hpp \
grammar_mutable.hpp \
grammar_simple.hpp \
grammar_static.hpp \
grammar_unknown.hpp \
graphviz.hpp \
head_finder.hpp \
hypergraph.hpp \
inside_outside.hpp \
intersect.hpp \
kbest.hpp \
lattice.hpp \
matcher.hpp \
model.hpp \
ngram.hpp \
ngram_index.hpp \
ngram_count_set.hpp \
operation.hpp \
operation_set.hpp \
optimize.hpp \
parameter.hpp \
parse.hpp \
parse_agenda.hpp \
parse_cky.hpp \
parse_coarse.hpp \
parse_tree.hpp \
parse_tree_cky.hpp \
permute.hpp \
prune.hpp \
prune_beam.hpp \
prune_density.hpp \
prune_edge.hpp \
prune_kbest.hpp \
push_bos_eos.hpp \
quantizer.hpp \
remove_annotation.hpp \
remove_epsilon.hpp \
rule.hpp \
semiring.hpp \
sentence.hpp \
sentence_vector.hpp \
signature.hpp \
sort.hpp \
sort_tail.hpp \
sort_topologically.hpp \
span_node.hpp \
span_edge.hpp \
span_forest.hpp \
span_vector.hpp \
statistics.hpp \
stemmer.hpp \
tokenizer.hpp \
treebank.hpp \
symbol.hpp \
symbol_vector.hpp \
transducer.hpp \
tree_grammar.hpp \
tree_grammar_mutable.hpp \
tree_grammar_static.hpp \
tree_grammar_simple.hpp \
tree_rule.hpp \
tree_transducer.hpp \
unite.hpp \
viterbi.hpp \
vocab.hpp \
weight_vector.hpp

libcicada_la_SOURCES = \
alignment.cpp \
attribute.cpp \
attribute_vector.cpp \
cluster.cpp \
feature.cpp \
feature_function.cpp \
format.cpp \
global_lexicon.cpp \
grammar_format.cpp \
grammar_mutable.cpp \
grammar_static.cpp \
grammar_unknown.cpp \
graphviz.cpp \
head_finder.cpp \
hypergraph.cpp \
lattice.cpp \
matcher.cpp \
model.cpp \
ngram.cpp \
ngram_index.cpp \
operation.cpp \
operation_set.cpp \
parameter.cpp \
rule.cpp \
sentence.cpp \
sentence_vector.cpp \
signature.cpp \
span_vector.cpp \
statistics.cpp \
stemmer.cpp \
tokenizer.cpp \
treebank.cpp \
symbol.cpp \
symbol_vector.cpp \
transducer.cpp \
tree_grammar_mutable.cpp \
tree_grammar_static.cpp \
tree_rule.cpp \
tree_transducer.cpp \
vocab.cpp

libcicada_la_LDFLAGS = \
	$(PTHREAD_FLAGS) \
	$(BOOST_THREAD_LDFLAGS) \
	$(BOOST_FILESYSTEM_LDFLAGS) \
	$(BOOST_IOSTREAMS_LDFLAGS) \
	$(BOOST_REGEX_LDFLAGS) \
	$(BOOST_SYSTEM_LDFLAGS) \
	-version-info $(CICADA_LTVERSION)

libcicada_la_LIBADD = \
	libcicada-eval.la \
	libcicada-format.la \
	libcicada-head.la \
	libcicada-feature.la \
	libcicada-matcher.la \
	libcicada-optimize.la \
	libcicada-operation.la \
	libcicada-semiring.la \
	libcicada-signature.la \
	libcicada-stemmer.la \
	\
	$(LIBUTILS) \
	$(BOOST_THREAD_LIBS) \
	$(BOOST_FILESYSTEM_LIBS) \
	$(BOOST_IOSTREAMS_LIBS) \
	$(BOOST_REGEX_LIBS) \
	$(BOOST_SYSTEM_LIBS) \
	$(ICU_LDFLAGS)

noinst_PROGRAMS = \
alignment_main \
attribute_vector_main \
cluster_main \
eval_main \
feature_vector_main \
format_main \
grammar_mutable_main \
grammar_static_main \
grammar_unknown_main \
hypergraph_main \
lattice_main \
matcher_main \
parameter_main \
rule_main \
sentence_main \
sentence_vector_main \
signature_main \
span_vector_main \
stemmer_main \
symbol_main \
tokenizer_main \
transducer_main \
tree_grammar_mutable_main \
tree_grammar_static_main \
tree_rule_main

alignment_main_SOURCES = alignment_main.cpp
alignment_main_LDADD = libcicada.la

attribute_vector_main_SOURCES = attribute_vector_main.cpp
attribute_vector_main_LDADD = libcicada.la

cluster_main_SOURCES = cluster_main.cpp
cluster_main_LDADD = libcicada.la

eval_main_SOURCES = eval_main.cpp
eval_main_LDADD = libcicada.la

feature_vector_main_SOURCES = feature_vector_main.cpp symbol.cpp vocab.cpp feature.cpp
feature_vector_main_CPPFLAGS = $(AM_CPPFLAGS)
feature_vector_main_LDADD = $(LIBUTILS)

format_main_SOURCES = format_main.cpp
format_main_LDADD = libcicada.la

grammar_mutable_main_SOURCES = grammar_mutable_main.cpp
grammar_mutable_main_LDADD = libcicada.la

grammar_static_main_SOURCES = grammar_static_main.cpp
grammar_static_main_LDADD = libcicada.la

grammar_unknown_main_SOURCES = grammar_unknown_main.cpp
grammar_unknown_main_LDADD = libcicada.la

hypergraph_main_SOURCES = hypergraph_main.cpp
hypergraph_main_LDADD = libcicada.la

lattice_main_SOURCES = lattice_main.cpp
lattice_main_LDADD = libcicada.la

matcher_main_SOURCES = matcher_main.cpp
matcher_main_LDADD = libcicada.la

parameter_main_SOURCES = parameter_main.cpp
parameter_main_LDADD = libcicada.la

rule_main_SOURCES = rule_main.cpp
rule_main_LDADD = libcicada.la

sentence_main_SOURCES = sentence_main.cpp
sentence_main_LDADD = libcicada.la

sentence_vector_main_SOURCES = sentence_vector_main.cpp
sentence_vector_main_LDADD = libcicada.la

signature_main_SOURCES = signature_main.cpp
signature_main_LDADD = libcicada.la

span_vector_main_SOURCES = span_vector_main.cpp
span_vector_main_LDADD = libcicada.la

stemmer_main_SOURCES = stemmer_main.cpp
stemmer_main_LDADD = libcicada.la

symbol_main_SOURCES = symbol_main.cpp
symbol_main_LDADD = libcicada.la

tokenizer_main_SOURCES = tokenizer_main.cpp
tokenizer_main_LDADD = libcicada.la

transducer_main_SOURCES = transducer_main.cpp
transducer_main_LDADD = libcicada.la

tree_grammar_mutable_main_SOURCES = tree_grammar_mutable_main.cpp
tree_grammar_mutable_main_LDADD = libcicada.la

tree_grammar_static_main_SOURCES = tree_grammar_static_main.cpp
tree_grammar_static_main_LDADD = libcicada.la

tree_rule_main_SOURCES = tree_rule_main.cpp
tree_rule_main_LDADD = libcicada.la
