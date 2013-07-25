#!/bin/sh

cicada=../../..

exec $cicada/scripts/cicada-index.py \
	--root-dir . \
	--model-dir . \
	\
	--ghkm \
	--sigtest-inclusive \
	\
	--feature-lexicon \
	--feature-root \
	--feature-internal \
	--feature-height \
	--threads 4

