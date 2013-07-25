#!/bin/sh

cicada=../../..

exec $cicada/scripts/cicada-index.py \
	--root-dir . \
	--model-dir . \
	\
	--scfg \
	--sigtest-inclusive \
	\
	--feature-lexicon \
	\
	--threads 4 

