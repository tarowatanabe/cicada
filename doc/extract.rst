Phrase Pair/Synchronous CFG/Synchronous TSG Extraction
======================================================

Phrase table or synchronous grammars, either synchronous context free
grammar (SCFG) or synchronous tree substitution grammar (STSG), are
learned by running the `cicada-extract.py` script.

.. code:: bash

  cicada-extract.py \
	  --f <source data> \
	  --e <target data> \
	  --ff <source hypergraph data> \
	  --fe <target hypergraph data> \
	  --a <word alignment> \
	  --{phrase, scfg, ghkm, tree}

``--f`` and ``--e`` options specify the source and target data, and
``--ff`` and ``--fe`` are their corresponding hypergraph data. The
format for hypergraph is documented in `doc/hypergraph.rst`.
Based on the option, you can extract different models:

--phrase  Extract phrase pairs.
--scfg    Extract synchronous-CFG.
--ghkm    Extract tree-to-string or string-to-tree model of
          synchronous-TSG. The model is determined whether we have
	  ``--ff`` option (tree-to-string) or ``--fe`` option
	  (string-to-tree).
--tree    Extract tree-to-tree model of synchronous-TSG. Both of
          ``--ff`` and ``--fe`` should be specified.

Details
-------

Model extraction is performed in three steps:

1. Lexical probabilities are estimated from the word aligned bilingual
   data for feature function computation. You can perform smoothing by

   - Variational Bayes estiamtes (``--lexicon-variational``, which is recommended)
   - L0 regularization (``lexicon-l0`` which is slower)

   The lexical probabilities `lex.f2n.gz` and `lex.n2f.gz` are dump in
   the directory `model` specified by the ``--lexical-dir`` options.

2. Extract phrase pairs or synchronous rules.
   Extracted counts are put on the subdirectory,
   `model/[grammar]-counts` where ``[grammar]`` is the model name,
   such as ``scfg`` or ``ghkm``. The `model` directory prefix can be
   modified by the ``--counts-dir`` option.
   
3. Merge extracted phrase pairs and synchronous rules.
   merged counts are put on the subdirectory,
   `model/[grammar]-score` where ``[grammar]`` is the model name,
   such as ``scfg`` or ``ghkm``. The `model` directory prefix can be
   modified by the ``--score-dir`` option.

Note that the merged model is not directly usable as a model, and the
extracted grammar should be further interpreted as a model. For
details, see `doc/indexing.hml`.
The format for the extracted and merged counts are summarized in
`man/cicada-extract.rst`.

The second and the third steps are efficiently computed by the
MapReduce framework using either pthreads or MPI, controlled by:

--threads        # of threads
--mpi            # of MPI jobs
--mpi-host       comma delimited list of hosts for MPI jobs
--mpi-host-file  host file for use with MPI jobs
--pbs            Run under pBS
--pbs-queue      PBS queue name

Extraction and final merging requires temporary disk space which can
be specified by 

--temporary-dir  temporary directory

Options
-------

There exists a couple of options which affect the phrase-table/grammar extractions.

SCFG
````

The SCFG in Moses can be simulated by the following options:

::

	--max-span-source 15
	--max-span-target 15
	--min-hole-source 2
	--min-hole-target 1
	--max-length 5
	--max-fertility 0
	--exhaustive 

The Hiero-style extraction from Chiang (2007) will be:

::

	--max-span-source 10
	--max-span-target 0
	--min-hole-source 2
	--min-hole-target 1
	--max-length 5
	--max-fertility 0
	--constrained

The ``--max-span-source`` and ``--max-span-target`` constrained the
length of the phrase pairs from which we will punch "holes" to
instantiate non-terminals. 0 implies no limits. ``--max-length`` is
the maximum number of terminals in the extracted synchronous
rules. Similarly, ``--max-fertility`` limits the maximum length ratio.
By default, we will extract all the phrase pairs, but uses minimum
sub-phrases to punch holes. This is adjusted by ``--exhaustive``
which uses all the sub-phrases as holes. The ``--constrained`` option
limits so that only minimum phrases are extracted from which we will
use sub-phrases to punch holes.

GHKM and Tree
`````````````

Here is a default setting in cicada:

::

	--max-nodes 7
	--max-height 4
	--max-compose 4

in which the maximum number of nodes in an elementary tree is 7
(``--max-nodes``), maximum height is 4 (``--max-height``) and the
number of minimum rules to compose a larger rule is limited to 4
(``--max-compose``).
There exist other options which greatly affect the extracted grammar:

--constrained  By default, all the minimum rules are extracted from
               the bilingual data, even if the minimum rules do not
	       satisfy other constraints, such as maximum number of
	       nodes or maximum height. This option force the minimum
	       rules to satisfy such constraints.
--exhaustive   By default, the non-aligned words are attached to the
               hyperedges, which are closer to the goal node. This
	       option try to attach non-aligned words exhaustively to
	       any possible hyperedges, and extract spuriously many
	       rules.
--project      In tree-to-string or string-to-tree models, the string
               side does not contain any linguistic labels. This
	       option project the non-terminal symbols from the
	       syntactic tree so that the string side rule contains
	       syntactic information. This is usually RECOMMENDED.
