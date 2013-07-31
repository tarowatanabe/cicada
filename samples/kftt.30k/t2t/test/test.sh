#!/bin/sh

cicada=../../../..
weights=

if test "$weights" = ""; then
  echo "where is your weights? (for instance ../tune/learn.10.weights)"
  exit 1
fi

### generate config file

$cicada/progs/cicada_filter_config \
  --input ../tune/cicada.config \
  --output cicada.config \
  --weights "weights=$weights" \
  --kbest 1 \
  --file "file=-"

### perform translation
$cicada/progs/cicada --config cicada.config --threads 4 --debug --input ../../t2s/data/dev.forest.ja.gz > dev.ja-en

### evaluation
$cicada/progs/cicada_eval --tstset dev.ja-en --refset ../../data/dev.en.ref


