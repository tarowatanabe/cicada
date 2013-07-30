Training
========

This is a step-by-step guide for training in cicada. We use a small
Japanese-English bilingual data distributed with the toolkit at
`samples/kftt.30k`. The data comes from the 30K sentences from
Kyoto Free Translation Task (http://www.phontron.com/kftt/), which was
derived from http://alaginrc.nict.go.jp/WikiCorpus/index_E.html.

Preprocessing
-------------

In this example, the Japanese data in `samples/kftt.30k/data/{train,tune,dev}.ja`
are segmented by Mecab (https://code.google.com/p/mecab/). Their
English parts are tokenized by a tokenizer in Moses (http://statmt.org/moses/)
with additional postprocessing to match the Penntreebank standard.

Word alignment
--------------

The first step in training is to align words in bilingual data. Cicada
can perform word alignment without relying on other tools.

.. code:: bash

  cd samples/kftt.30k/aignment
  ../../../scripts/cicada-alignment.py \
	--f ../data/train.ja.bz2 \
	--e ../data/train.en.bz2 \
	--symmetric \
	--variational \
	--posterior \
	--threads 4

In this example, by default, we run 5 of IBM Model 1, 5 of HMM, 5 of
IBM Model 4 exploiting 4 threads (**--threads 4**). During the training,
posteriors are constrained (**--posterior**) so that bilingual words
are symmetized (**--symmetric**), and the parameters are smoothe by
a naive variational Bayes estiamte (**--variational**). By default,
the word alignment data is found at `samples/kftt.30k/aignment/model/aligned.posterior-itg`
which is word alignment constrained by ITG. For details, see
`doc/alignment.rst`.

Translation Model
-----------------

The second step of training is an estimation of translation model. In
this example, we will construct a synchronous-CFG, SCFG, also known as
Hiero grammar. If you want to train other models, such as
tree-to-string or string-to-tree models, see `doc/training-stsg.rst`.

The training model is constructed in two steps, firstly, collecting
synchronous rules, secondly, estimating features for each synchronous
rule.

Extraction
``````````
.. code:: bash

  cd samples/kftt.30k/scfg/model
  ../../../../scripts/cicada-extract.py \
	--f ../../data/train.ja.bz2 \
	--e ../../data/train.en.bz2 \
	--a ../../alignment/model/aligned.posterior-itg \
	--model-dir . \
	--lexicon-variational \
	--scfg \
	--exhaustive \
	--threads 4

In this example, the SCFG model (**--scfg**) is output at the current
directory (**--model-dir .**), and lexicon model is also estimated
with a naive variational Bayes estimate (**--lexicon-variational**).
The rules are exhaustively extracted (**--exhaustive**) by usign all
the phrase pairs as holes. Extracted counts are stored at
`samples/kftt.30k/scfg/model/scfg-counts`.
Since we exploited 4 threads (**--threads 4**), 4 model files
``[0-3].gz`` are created at `samples/kftt.30k/scfg/model/scfg-score`.
For details, see `doc/extract.rst`.

Features
````````

The models created at `samples/kftt.30k/scfg/model/scfg-score`
contain only count information, and should be interpreted as features.

.. code:: bash

  ../../../../scripts/cicada-index.py \
	--model-dir . \
	--scfg \
	--sigtest-inclusive \
	--feature-lexicon \
	--threads 4 

This is an example to generate 4 models ``[0-3].gz`` and corresponding
indexed binary models ``[0-3].bin`` at `samples/kftt.30k/scfg/model/scfg-index`
using 4 threads (**--threads 4**). We exclude synchronous rules which
are less significant measured by Fisher's exact test, but preserves rare
rules (**--sigtest-inclusive**). By default, we use 2 features, which
are :math:`\log p(source | target)` and :math:`\log p(target | source)`
estimated by the relative frequencies. **--feature-lexicon** option
includes lexical weight features in two directions.

Language Model
--------------

The ngram language model is one of the important models in machine
translation.  There exists many alternative options to train ngram
language models, but in this example, we use `expgram`, an ngram
toolkit available from http://www2.nict.go.jp/univ-com/multi_trans/expgram.

.. code:: bash

  cd samples/kftt.30k/ngram
  [directory-for-expgram]/progs/expgram_counts_extract \
    --corpus ../data/train.en.bz2 \
    --output ngram.5.en.counts \
    --order 5 \
    --threads 4
  [directory-for-expgram]/progs/expgram_counts_estimate \
    --ngram ngram.5.en.counts \
    --output ngram.5.en.lm \
    --shard 4

Here, we use 4 threads to estimate an ngrma language model by, first,
collecting counts (**expgram_counts_extract**), then, by estimating
the model (**expgram_counts_estimate**). An alternative is to use a
script included in the expgram:

.. code:: bash

  cd samples/kftt.30k/ngram
  [directory-for-expgram]/scripts/expgram.py \
    --corpus ../data/train.en \
    --output ngram.5.en \
    --threads 4

Tuning
------

Now, we are ready to perform translation, but it is better to tune
the parameters to determine the combination weights for features, such
as translation models and language models.

Configuration file
``````````````````

First, we need to create a configuration file to run decoder.

.. code:: bash

  cd samples/kftt.30k/tune
  ../../../../scripts/cicada-config.py \
	--grammar ../model/scfg-index \
	--max-span 15 \
	--straight \
	--insertion \
	--feature-ngram ../../ngram/ngram.5.en.lm \
	--scfg \
	--beam 200 > cicada.config

In this example, we use the grammar in `..model/scfg-index` with
maximum span set to 15 (**--max-span 15**). As glue rules, we employ
monotone rule (**--straight**) and use insertion grammar to copy the
input string into output string (**--insertion**). We use additional
ngram language model feature (**--feature-ngram**) in the
model. Translation is carried out by SCFG decoding (**--scfg**) with
beam size of 200 (**--beam 200**) for cube pruning.

The configuration file consists of 3 parts, grammars, features and
operations. For details, see `doc/grammar.rst`, `doc/features.rst` and
`doc/operation.rst`. Actually, the configuration file is a template so
that we can instantiate parameters ``${weights}`` and output files
``${file}`` during tuning or when testing.

Reference translations
``````````````````````

During tuning, we need reference translations a set of high quality
translations for each input sentence. In cicada, multiple refernce
translations are summarized in a single file as follows:

::

   0 ||| first reference
   0 ||| second refernece
   1 ||| first reference for the second input
   1 ||| second reference for the second input

The format is very simple and we also provide a program to generate
the file from multiple translations:

.. code:: bash

  cd samples/kftt.30k/data
  ../../../progs/cicada_filter_refset tune.en --output tune.en.ref
  ../../../progs/cicada_filter_refset dev.en --output dev.en.ref

Tune parameters
```````````````

Now, we are ready to perfom tuning:

.. code:: bash

  cd samples/kftt.30k/scfg/tune
  ../../../../scripts/cicada-learn.py \
	--srcset ../../data/tune.ja \
	--refset ../../data/tune.en.ref \
	--config cicada.config \
	--kbest 1000 \
	--threads 4

We use `tune.ja` as a source set (**--srcset**) and `tune.en.ref` as
its reference translations (**--refset**) with `cicada.config` as a
configuratin template (**--config**). The training is performed by
k-best merging batch style learning with 1,000 best transaltions
generated in each round (**--kbest**). By default, training objective
is xBLEU, which is superior to other objectives, like pair-wise
ranking (PRO) or direct error minimization (MERT).

By default, training is performed 10 iterations, and generates several
files ``learn.<iteration>.*``. The tuned parameters have suffix of
``.weights``.

Testing
-------

After tuning, we can perform actual translation for test data. First,
we will generate a configuration file from the template:

.. code:: bash

  cd samples/kftt.30k/scfg/test
  ../../../../progs/cicada_filter_config \
    --input ../tune/cicada.config \
    --output cicada.config \
    --weights "weights=../tune/learn.10.weights" \
    --kbest 1 \
    --file "file=-"

In this example, we use exactly the same template employed for tuning
(**--input**), and use the parameters learned after the last iteartion
(**--weights**). The single-best (**--kbest 1**) is output to stdout
(**--file**).

Decoding
````````

Then, we translate the development data `dev.ja` and output as `dev.ja-en`:

.. code:: bash

  ../../../../progs/cicada \
	  --config cicada.config \
	  --threads 4  < ../../data/dev.ja > dev.ja-en

Evalluation
```````````

Finally, we will evaluate the translation:

.. code:: bash

  ../../../../progs/cicada_eval \
	  --tstset dev.ja-en \
	  --refset ../../data/dev.en.ref

which computes BLEU score, by default. For details of evaluation, see
`doc/eval.rst`.
