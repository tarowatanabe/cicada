#!/bin/sh

cicada=../..

exec $cicada/progs/cicada \
      --input $cicada/samples/scfg/input.txt \
      --grammar $cicada/samples/scfg/grammar.bin \
      --grammar "glue:straight=true,inverted=false,non-terminal=[x],goal=[s]" \
      --grammar "insertion:non-terminal=[x]" \
      --feature-function "ngram:file=$cicada/samples/scfg/ngram.bin" \
      --feature-function word-penalty \
      --feature-function rule-penalty \
      --operation compose-cky \
      --operation apply:prune=true,size=100,weights=$cicada/samples/scfg/weights \
      --operation output:file=-,kbest=10,weights=$cicada/samples/scfg/weights
