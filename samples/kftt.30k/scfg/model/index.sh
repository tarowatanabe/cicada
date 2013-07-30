#!/bin/sh

cicada=../../../..

exec $cicada/scripts/cicada-index.py \
	--model-dir . \
	\
	--scfg \
	--sigtest-inclusive \
	\
	--feature-lexicon \
	\
	--threads 4 

