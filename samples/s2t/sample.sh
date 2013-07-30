#!/bin/sh

cicada=../..

##
##
##    gramamr:
##      tree-grammar at $cicada/samples/s2t/grammar.bin
##      insertion grammar using [x] as non-terminal
##      The target side goal symbol is [ROOT]
##    features:
##      ngram language model
##      # of words
##      # of rules
##      # of glue penalty which penalize the # of SCFG rule to be included in the forest
##    operations:
##      tree-cky-composition
##      apply features
##      output kbests 

exec $cicada/progs/cicada \
      --input $cicada/samples/scfg/input.txt \
      --goal '[ROOT]' \
      --tree-grammar $cicada/samples/s2t/grammar.bin \
      --grammar "insertion:non-terminal=[x]" \
      --feature-function "ngram:file=$cicada/samples/scfg/ngram.bin" \
      --feature-function word-penalty \
      --feature-function rule-penalty \
      --feature-function glue-tree-penalty \
      --operation compose-tree-cky \
      --operation apply:prune=true,size=100,weights=$cicada/samples/s2t/weights \
      --operation output:file=-,kbest=10,weights=$cicada/samples/s2t/weights
