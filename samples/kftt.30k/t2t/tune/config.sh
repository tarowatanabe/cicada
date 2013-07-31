#!/bin/sh

cicada=../../../..

exec $cicada/scripts/cicada-config.py \
	--tree-grammar ../model/tree-index \
   	--goal '[ROOT]' \
        --glue '[x]' \
	--fallback-glue \
	--feature-ngram ../../ngram/ngram.5.en.lm \
	--tree \
	--beam 1024 > cicada.config
	
