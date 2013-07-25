#!/bin/sh

cicada=../../..

exec $cicada/scripts/cicada-extract.py \
	--f ../../data/train.ja \
	--e ../../data/train.en \
	--a ../../alignment/model/aligned.posterior-itg \
	--ff ../data/train.forest.ja.gz \
	\
	--root-dir . \
	--model-dir . \
	\
	--lexicon-variational \
	\
	--ghkm \
	--constrained \
	\
	--threads 4 \
	--max-malloc 16

