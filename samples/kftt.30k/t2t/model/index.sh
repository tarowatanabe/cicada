#!/bin/sh

cicada=../../../..

exec $cicada/scripts/cicada-index.py \
	--model-dir . \
	\
	--tree \
	--sigtest-inclusive \
	\
	--feature-lexicon \
	--feature-root \
	--feature-internal \
	--feature-height \
	--threads 4

