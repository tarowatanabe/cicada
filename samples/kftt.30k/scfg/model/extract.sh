#!/bin/sh

cicada=../../..

exec $cicada/scripts/cicada-extract.py \
	--f ../../data/train.ja.bz2 \
	--e ../../data/train.en.bz2 \
	--a ../../alignment/model/aligned.posterior-itg \
	\
	--root-dir . \
	--model-dir . \
	\
	--lexicon-variational \
	\
	--scfg \
	--exhaustive \
	\
	--threads 4 \
	--max-malloc 32

