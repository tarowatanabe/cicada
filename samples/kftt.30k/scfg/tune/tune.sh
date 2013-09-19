#!/bin/sh

cicada=../../../..

exec $cicada/scripts/cicada-learn.py \
	--srcset ../../data/tune.ja \
	--refset ../../data/tune.en.ref \
	--config cicada.config \
	--kbest 1000 \
	--regularize-l2 1e-5 \
        --merge \
	--threads 4
