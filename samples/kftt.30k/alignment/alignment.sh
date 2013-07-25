#!/bin/sh

cicada=../../..

## We perform symmetized posterior constrained training, and
## perform smoothing by naive Bayes.

exec ${cicada}/scripts/cicada-alignment.py \
	--f ../data/train.ja.bz2 \
	--e ../data/train.en.bz2 \
	--symmetric \
	--variational \
	--posterior \
	--threads 4
