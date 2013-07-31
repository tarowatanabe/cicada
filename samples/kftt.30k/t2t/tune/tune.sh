#!/bin/sh

cicada=/home/t_watana/research/cicada

exec $cicada/scripts/cicada-learn.py \
	--srcset ../../t2s/data/tune.forest.ja.gz \
	--refset ../../data/tune.en.ref \
	--config cicada.config \
	--kbest 1000 \
	--threads 16
