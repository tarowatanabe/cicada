#!/bin/sh

cicada=../../../..

exec $cicada/scripts/cicada-config.py \
	--tree-grammar ../model/ghkm-index \
	--goal '[x]' \
	--glue '[x]' \
	--fallback-glue \
	--feature-ngram ../../ngram/ngram.5.en.lm \
	--ghkm \
	--beam 1024 > cicada.config
	
