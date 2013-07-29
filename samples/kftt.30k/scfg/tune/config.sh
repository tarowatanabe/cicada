#!/bin/sh

cicada=../../../..

exec $cicada/scripts/cicada-config.py \
	--grammar ../model/scfg-index \
	--max-span 15 \
	--straight \
	--insertion \
	--feature-ngram ../../ngram/ngram.5.en.lm \
	--scfg \
	--beam 200 > cicada.config
	
