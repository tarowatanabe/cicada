#!/bin/sh

cicada=../../../..

exec $cicada/scripts/cicada-extract.py \
	--f ../../data/train.ja \
	--e ../../data/train.en \
	--a ../../alignment/model/aligned.posterior-itg \
	--ff ../../t2s/data/train.forest.ja.gz \
	--fe ../../s2t/data/train.tree.en.gz \
	\
	--model-dir . \
	\
	--lexicon-variational \
	\
	--tree \
	--constrained \
        --exhaustive \
	\
	--threads 4

