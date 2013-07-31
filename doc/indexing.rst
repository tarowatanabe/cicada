Indexing for Phrase Pair/Synchronous CFG/Synchronous TSG
========================================================

After model extraction from bilingual data, we generate a synchronous
grammar with features and also index into a binary format for
efficient query. The procedure is performed by the `cicada-index.py`
script.

.. code:: bash

  cicada-index.py --{phrase, scfg, ghkm, tree}

The option corresponds to the extracted model by the
`cicada-extract.py`. When you translate using a string-to-tree model,
you need an additional flag ``--cky`` to perform indexing suitable for
the CKY algorithm.
All the counts in `model/[grammar]-counts` and `model/[grammar]-score`
are reinterpreted and dump into `model/[grammar]-index` in both of
plain gzip file and a binary format.

Pruning
-------

We support various grammar pruning method. One of the best methods is
``--sigtest-inclusive`` or ``--sigtest-exclusive`` which has been
proved useful especially in extremely large bilingual data.

--kbest    Retain only kbest of synchronous rules measured by frequencies.
--sigtest  Retain only statistically significant synchronous rules
           according to Fisher's exact test.
	   ``--sigtest-inclusive`` keep rare rules which occurs only
	   once, but ``--sigtest-exclusive`` discards such rules.
--cutoff   Count cutoff.
--threshold  Probability threshold by :math:`p(source|target)` and :math:`p(target|source)`

Features
--------

All the features are indexed as a float value, but you can perform
8-bit quantization for compact storage via ``--quantize`` option.
By default, only two features are encoded in the model:
:math:`\log p(source|target)` and :math:`\log p(target|source)`.
Additional features are included by the following options:

--feature-root        generative probability
--feature-fisher      Fisher's exact test
--feature-type        observation probability
--feature-singleton   singleton features
--feature-cross       cross features
--feature-unaligned   unaligned features
--feature-internal    internal features
--feature-height      height features
--feature-lexicon     compute Model1 features
--feature-model1      compute Model1 features
--feature-noisy-or    compute noisy-or features
--feature-insertion-deletion  Insertion deletion features

