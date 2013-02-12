Feature functions in cicada

You can assign features to each rule by indicating feature names for
each rule:

[x] ||| [x,1] good ||| [x,1] bad ||| feature1=0.5 another-feature5=0.6

or, when specifying by --grammar option, you can include additional parameters,

<grammar-file>:feature0=a-name,feature1=another-name, etc.

which will assign non-named features to pre-specified names. By default, we will
assign, rule-table-0, rule-table-1, etc.

For "static grammar", such as grammar learned by moses, joshua, nicttm etc.,
you cannot assign arbitrary different features, since the number of features
should be equal. However, you can assign feature-names, as in non-static mutable
grammar:
	--grammar-static <grammar-file>:feature0=name,feature1=name1

Other dynamically assigned features are applied via "apply" operaion:

antecedent: sparse features assigned two-level non-terminals, delimited terminals and span-size
bleu: bleu score
blue-linear: linear bleu score
bleu-expected: expected bleu score (You need to compute expected counts before this feature application)
bleu-multi: multiple bleu score feature
neighbours: sparse features assigned for non-terminal, its neighboring words and span-size
ngram: ngram langauge model feature. You can assign multiple LM by defining name=<feature-name> option.
ngram-tree: sparse featuers assigned for non-terminal path for boundaries.
parent:  rule (itself) and rule with its parent category
permute: permutation feature (we assume tree is permuted via "permute" operation)
rule-shape: rule-shape feature, which book keep the permuted non-terminal + terminal
sgml-tag: SGML tag feature
span: span feature
sparse-lexicon: sparse lexicon feature
sparse-ngram: sparse ngram feature
variational: features for variational decoding
vocabulary: in-vocabulary feature
arity-penalty: arity penalty feature
glue-tree-penalty: glue-tree penalty
non-latin-penalty: non-latin penalty
word-penalty: simply count # of words
rule-penalty: simply count # of rules

HOW TO DEVELOP YOUR FEATURE FUNCTIONS?

The feature function classes are located at cicada/feature
Basically, you have to inherit FeatureFunction class in cicada/feature_function.hpp
and fill-in virtual function. See the comments.

Then, you want to add initializer in cicada/feature_function.cpp so that
FeatureFunction::create(string) can generate correct feature funtion.
