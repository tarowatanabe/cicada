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

###
### for tuning and development data only
###
### Currently, string-to-tree model is limited in that we have no hard limit in
### translation. Thus, the decoding is extremely slow and consumes large memory.
### Here, as a sample, we will limi the sentence length to 20.
### This will be resolved in the near future by enforcing span-limit as in SCFG.

for data in tune dev; do
  ./length-limit.py \
    --source ../../data/$data.ja \
    --target ../../data/$data.en.ref \
    --output-source $data.ja \
    --output-target $data.en.ref
done
