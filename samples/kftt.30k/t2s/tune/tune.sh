#!/bin/sh

cicada=../../..

exec $cicada/scripts/cicada-learn.py \
	--devset ../data/tune.forest.ja.gz \
	--refset ../../data/tune.en.ref \
	--config cicada.config \
	--kbest 1000 \
	--threads 4
