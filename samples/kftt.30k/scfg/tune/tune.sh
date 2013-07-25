#!/bin/sh

cicada=../../..

exec $cicada/scripts/cicada-learn.py \
	--devset ../../data/tune.ja \
	--refset ../../data/tune.en.ref \
	--config cicada.config \
	--kbest 1000 \
	--threads 16
