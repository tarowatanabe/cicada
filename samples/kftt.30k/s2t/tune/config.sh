#!/bin/sh

cicada=../../../..

exec $cicada/scripts/cicada-config.py \
	--tree-grammar ../model/ghkm-index \
	--max-span 10 \
	--goal '[ROOT]' \
        --tree-straight \
	--insertion \
	--feature-ngram ../../ngram/ngram.5.en.lm \
	--tree-cky \
	--beam 100 > cicada.config
	
