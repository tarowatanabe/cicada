#!/bin/sh

cicada=../../../..

exec $cicada/scripts/cicada-config.py \
	--tree-grammar ../model/ghkm-index \
	--max-span 20 \
	--goal '[ROOT]' \
	--insertion \
	--feature-ngram ../../ngram/ngram.5.en.lm \
	--tree-cky \
	--beam 200 > cicada.config
	
