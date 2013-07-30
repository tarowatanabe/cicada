#!/bin/sh

cicada=../..

##
## 1. pentreebank is converted into hypergraph
##
## 2. Then, translate
##    gramamr:
##      tree-grammar at $cicada/samples/t2s/grammar.bin
##      fallback grammar
##      The target side goal symbol is [ROOT]
##    features:
##      ngram language model
##      # of words
##      # of rules
##    operations:
##      cyk-binarization
##      tree-composition
##      apply features
##      output kbests 

$cicada/progs/cicada_filter_penntreebank \
      --input $cicada/samples/t2s/input.txt \
      --normalize \
| \
$cicada/progs/cicada \
      --input - \
      --input-forest \
      --goal '[ROOT]' \
      --tree-grammar $cicada/samples/t2s/grammar.bin \
      --tree-grammar fallback \
      --feature-function "ngram:file=$cicada/samples/scfg/ngram.bin" \
      --feature-function word-penalty \
      --feature-function rule-penalty \
      --operation binarize:direction=cyk,order=1 \
      --operation compose-tree \
      --operation apply:prune=true,size=1000,weights=$cicada/samples/t2s/weights \
      --operation output:file=-,kbest=10,weights=$cicada/samples/t2s/weights
