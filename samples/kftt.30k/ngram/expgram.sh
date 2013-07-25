#!/bin/sh

expgram={where did you installed expgram?}

### Following is a quick example of LM estimation on a small data set.

$expgram/progs/expgram_counts_extract --corpus ../data/train.en --output ngram.5.en.counts --order 5 --threads 2 --debug

$expgram/progs/expgram_counts_estimate --ngram ngram.5.en.counts --output ngram.5.en.lm --shard 2 --debug

### A standard way is:
# $expgram/scripts/expgram.py --corpus ../data/train.en --output ngram.5.en --threads 2
