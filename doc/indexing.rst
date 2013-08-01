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
the CKY algorithm. Note that tree-to-{tree,string} models can be
indexed by ``--cky`` when composing with strings, not trees.
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

--feature-root        Generative probabilities: :math:`\log p(source|root(source))` and :math:`\log p(target|root(target))`.
--feature-fisher      Fisher's exact test
--feature-type        observation probability
--feature-singleton   Singleton features: true if this is pair occurs 
                      only once in bilingual data.
--feature-cross       Cross features: the number of crossing in word alignment.
--feature-unaligned   Unaligned features: the number of words which
                      are not aligned with each other.
--feature-internal    Internal features (STSG only). The number of
                      internal nodes in an elementary tree.
--feature-height      Height features: the height of elementary trees.
--feature-lexicon     Lexical features
--feature-model1      Model1 features
--feature-noisy-or    Noisy-or features
--feature-insertion-deletion  Insertion deletion features

Parallel Indexing
-----------------

Similar to model extraction procedure, we support parallel indexing
using either pthreads or MPI, controlled by:

--threads        # of threads
--mpi            # of MPI jobs
--mpi-host       comma delimited list of hosts for MPI jobs
--mpi-host-file  host file for use with MPI jobs
--pbs            Run under PBS
--pbs-queue      PBS queue name

Indexing requires temporary space which can be specified by the
option:

--temporary-dir  temporary directory
