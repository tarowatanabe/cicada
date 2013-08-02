#!/bin/sh

stanford=
cicada=../../../..

if test "$stanford" = ""; then
  echo "where is your stanford parser?"
  exit 1
fi

# Here, we use stanford-parser to parse training data in English
# cicada_filter_penntreebank to transform into hypergraph.
#   Note that the stanford parser re-normalize words, like '(' into '-LRB-'.
#   we use --map <file name> to uncover the original tokenizatin.
#   The penntreebank has a very strange label, like ','.  Thus, we normalize constituency labels into COMMA etc.
# cicada to binarize in the left-heavy direction, and the labels memory only two contexts.

bzcat ../../data/train.en.bz2 | \
java \
    -mx12g \
    -cp $stanford/stanford-parser.jar:$stanford/stanford-parser-3.2.0-models.jar \
    -tLPP edu.stanford.nlp.parser.lexparser.EnglishTreebankParserParams \
    -tokenized -sentences newline \
    -escaper edu.stanford.nlp.process.PTBEscapingProcessor \
    -encoding UTF-8 \
    -maxLength 400 \
    -MAX_ITEMS 800000 \
    -outputFormat oneline \
    -outputFormatOptions includePunctuationDependencies \
    edu/stanford/nlp/models/lexparser/englishFactored.ser.gz \
    - | \
$cicada/progs/cicada_filter_penntreebank \
    --map ../../data/train.en.bz2 \
    --normalize | \
$cicada/progs/cicada \
    --input-forest \
    --threads 8 \
    --operation binarize:direction=left,order=2 \
    --operation output:no-id=true,file=train.tree.en.gz
