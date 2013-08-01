Feature functions
=================

Other dynamically assigned features are applied via "apply"
operations.

Features
--------

The list of available features can be found by ``cicada
--feature-function-list``, and summarized as follows:

antecedent: antecedent feature
	cluster=[word class file]
	stemmer=[stemmer spec]
	alignment=[true|false] alignment forest mode
	source-root=[true|false] source root mode (for tree composition)
	name=feature-name-prefix(default: antecedent)

bleu: BLEU
	order=<order>
	exact=[true|false] clipped ngram computation
	skip-sgml-tag=[true|false] skip sgml tags
	tokenizer=[tokenizer spec]
	name=feature-name(default: bleu)
	refset=reference set file

bleus: BLEUS
	order=<order>
	exact=[true|false] clipped ngram computation
	skip-sgml-tag=[true|false] skip sgml tags
	tokenizer=[tokenizer spec]
	name=feature-name(default: bleus)
	refset=reference set file

bleu-expected: expected-BLEU
	order=<order>

bleu-linear: linear corpus-BLEU
	order=<order>
	precision=<default 0.8>
	ratio=<default 0.6>
	tokenizer=[tokenizer spec]
	name=feature-name(default: bleu-linear)
	refset=reference set file
	skip-sgml-tag=[true|false] skip sgml tags

bleu-multi: multiple BLEU
	order=<order>
	exact=[true|false] clipped ngram computation
	skip-sgml-tag=[true|false] skip sgml tags
	tokenizer=[tokenizer spec]
	name=feature-name (default: bleu-multi:integer)
	size=# of BLEU features

deletion: deletion feature
	file=lexicon model file
	populate=[true|false] "populate" by pre-fetching
	skip-sgml-tag=[true|false] skip sgml tags
	unique-source=[true|source] unique source labels
	name=feature-name (default: deletion)
	threshold=threshold (default: 0.5)

dependency: dependency feature
	order=[order] (default: 2)

depeval: dependency evaluation feature
	skip-sgml-tag=[true|false] skip sgml tags
	tokenizer=[tokenizer spec]

distortion: phrase-based distortion

frontier-bigram: sparse frontier source side bigram
	source=[true|false] source side bigram (this is a default)
	target=[true|false] target side bigram (you can specify both)
	skip-sgml-tag=[true|false] skip sgml tags
	name=feature-name-prefix (default: frontier-bigram)

frontier-lexicon: sparse lexicon feature from frontiers
	cluster-source=[word class file] word-class for source side
	cluster-target=[word class file] word-class for target side
	stemmer-source=[stemmer spec] stemming for source side
	stemmer-target=[stemmer spec] stemming for target side
	skip-sgml-tag=[true|false] skip sgml tags
	name=feature-name-prefix (default: frontier-lexicon)

frontier-pair: sparse frontier pair features
	skip-sgml-tag=[true|false] skip sgml tags
	name=feature-name-prefix (default: frontier-pair)

frontier-shape: sparse frontier shape features
	skip-sgml-tag=[true|false] skip sgml tags
	name=feature-name-prefix (default: frontier-shape)

global-lexicon: global lexicon feature
	file=global lexicon file

insertion: insertion feature
	file=lexicon model file
	populate=[true|false] "populate" by pre-fetching
	skip-sgml-tag=[true|false] skip sgml tags
	unique-source=[true|source] unique source labels
	name=feature-name (default: insertion)
	threshold=threshold (default: 0.5)

kenlm: ngram LM from kenlm
	file=<file>
	populate=[true|false] "populate" by pre-fetching
	cluster=<word class>
	name=feature-name(default: ngram)
	no-bos-eos=[true|false] do not add bos/eos
	skip-sgml-tag=[true|false] skip sgml tags
	coarse-file=<file>   ngram for coarrse heuristic
	coarse-populate=[true|false] "populate" by pre-fetching
	coarse-cluster=<word class> word class for coarse heuristics

lexicalized-reordering: lexicalized reordering for phrase-based
	bidirectional=[true|false]
	monotonicity=[true|false]
	feature=attribute name mapping

lexicon: lexicon model feature based on P(target-sentence | source-sentence)
	file=lexicon model file
	populate=[true|false] "populate" by pre-fetching
	skip-sgml-tag=[true|false] skip sgml tags
	unique-source=[true|source] unique source labels
	name=feature-name-prefix (default: lexicon)

neighbours: neighbour words feature
	cluster=[word class file]
	stemmer=[stemmer spec]
	alignment=[true|false] alignment forest mode
	source-root=[true|false] source root mode (for tree composition)
	name=feature-name-prefix(default: neighbours)

ngram: ngram language model
	file=<file>
	populate=[true|false] "populate" by pre-fetching
	cluster=<word class>
	name=feature-name(default: ngram)
	approximate=[true|false] approximated upper-bound estimates
	no-bos-eos=[true|false] do not add bos/eos
	skip-sgml-tag=[true|false] skip sgml tags
	coarse-file=<file>   ngram for coarrse heuristic
	coarse-populate=[true|false] "populate" by pre-fetching
	coarse-cluster=<word class> word class for coarse heuristics
	coarse-approximate=[true|false] approximated upper-bound estimates

ngram-pyp: Pitman-Yor Process ngram language model
	file=<file>
	populate=[true|false] "populate" by pre-fetching
	order=<order>
	name=feature-name(default: ngram-pyp)
	no-bos-eos=[true|false] do not add bos/eos
	skip-sgml-tag=[true|false] skip sgml tags
	coarse-order=<order> ngram order for coarse heuristic
	coarse-file=<file>   ngram for coarrse heuristic
	coarse-populate=[true|false] "populate" by pre-fetching

ngram-tree: ngram tree feature
	cluster=[word class file]
	stemmer=[stemmer spec]
	alignment=[true|false] alignment forest mode
	source-root=[true|false] source root mode (for tree composition)
	name=feature-name-prefix(default: ngram-tree)

parent: parent feature
	cluster=[word class file]
	stemmer=[stemmer spec]
	exclude-terminal=[true|false] exclude terminal symbol
	alignment=[true|false] alignment forest mode
	source-root=[true|false] source root mode (for tree composition)
	name=feature-name-prefix(default: parent)

permute: permutation feature
	weights=weight file for collapsed feature
	collapse=[true|false] collapsed feature

sgml-tag: sgml-tag feature

span: lexical span feature

variational: variational feature for variational decoding
	order=<order>
	no-bos-eos=[true|false] do not add bos/eos

sparse-lexicon: sparse lexicon feature
	cluster-source=[word class file] word-class for source side
	cluster-target=[word class file] word-class for target side
	stemmer-source=[stemmer spec] stemming for source side
	stemmer-target=[stemmer spec] stemming for target side
	skip-sgml-tag=[true|false] skip sgml tags
	unique-source=[true|source] unique source labels
	name=feature-name-prefix (default: sparse-lexicon)
	lexicon=lexicon file
	lexicon-prefix=prefix lexicon file
	lexicon-suffix=suffix lexicon file
	pair=[true|false]   use of simple source/target pair (default: true)
	prefix=[true|false] use of source prefix
	suffix=[true|false] use of source suffix
	fertility=[true|false] lexical fertility feature

sparse-ngram: sparse ngram feature
	order=<order>
	no-bos-eos=[true|false] do not add bos/eos
	skip-sgml-tag=[true|false] skip sgml tags
	name=feature-name-prefix (default: sparse-ngram)
	cluster=[word class file] word-class for ngram
	stemmer=[stemmer spec] stemming for ngram

word-penalty: word penalty feature

rule-penalty: rule penalty feature

arity-penalty: rule arity penalty feature

glue-tree-penalty: glue tree penalty feature

non-latin-penalty: non-latin word penalty feature

internal-node-penalty: internal node penalty feature

rule-shape: rule shape feature

vocabulary: vocabulary feature
	file=[vocabulary file]
	oov=[true|false] oov penalty mode
	name=feature-name (default: vocabulary)

relative-position: relative alignment feature
	cluster=[word class file]
	stemmer=[stemmer spec]

path: path features
	cluster-source=[word class file]
	cluster-target=[word class file]
	stemmer-source=[stemmer spec]
	stemmer-target=[stemmer spec]

null-path: path involving epsilon

fertility-local: local fertility feature
	cluster=[word class file]
	stemmer=[stemmer spec]

target-bigram: target bigram feature
	cluster=[word class file]
	stemmer=[stemmer spec]

word-pair: word pair feature
	cluster-source=[word class file]
	cluster-target=[word class file]
	stemmer-source=[stemmer spec]
	stemmer-target=[stemmer spec]



Develop a new feature function
------------------------------

The feature function classes are located at `cicada/feature`
Basically, you have to inherit FeatureFunction class in `cicada/feature_function.hpp`
and fill-in virtual function. See the comments.

Then, you want to add initializer in cicada/feature_function.cpp so that
FeatureFunction::create(option) can generate correct feature functions.
