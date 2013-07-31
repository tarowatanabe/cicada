Training {tree,string}-to-{tree,string} Models
==============================================

This is a step-by-step guide for training syntactical models in
cicada. Following the `doc/training.rst`, we use the toy examples in
`samples/kftt.30k`.

Tree-to-string Model
--------------------

Preprocessing
`````````````

In tree-to-{string,tree} modeling, input is a parser tree, or
forest. As an example, we use CaboCha (http://code.google.com/p/cabocha/)
to dependency parse the Japanese side and transform the tree structure
into a constituency structure.

.. code:: bash

  cd samples/kftt.30k/t2s/data
  bzcat ../../data/train.ja.bz2 | \
  awk '{for (i=1;i<=NF;++i) {printf "%s\t*\n", $i } print "EOS";}' | \
  [directory-for-mecab]/mecab -p | \
  [directory-for-cabocha]/cabocha -f1 -I 1 | \
  ../../../../progs/cicada_filter_dependency \
	--cabocha \
	--func \
	--forest \
	--head | \
  ../../../../progs/cicada \
	--input-forest \
	--threads 4 \
	--operation binarize:direction=cyk,order=1 \
	--operation output:no-id=true,file=train.forest.ja.gz

In this example, Japanese is POS tagged by Mecab, dependency
parsed by CaboCha, then transformed into a hypergraph format
(``cicada_filter_dependency``). The dependency structure is
transformed into a constituency structure using function words as a
head (``--func``). The tree structure is exhaustively binarized by the
CYK binarization algorithm with order 1 (``binarize:direction=cyk,order=1``).
Other data, such as development and testing data (`tune.ja` and
`dev.ja`) are also preprocessed similarly.

Extraction
``````````

We use the word alignment extracted in `samples/kftt.30k/alignment`
and extract synchronous rules:

.. code:: bash

  cd samples/kftt.30k/t2s/model
  ../../../../scripts/cicada-extract.py \
	--f ../../data/train.ja.bz2 \
	--e ../../data/train.en.bz2 \
	--a ../../alignment/model/aligned.posterior-itg \
	--ff ../data/train.forest.ja.gz \
	--model-dir . \
	--lexicon-variational \
	--ghkm \
	--constrained \
	--threads 4

Major difference to the SCFG training is to specify the source side
forest (``--ff``) in addition to the Japanese and English bilingual
sentences (``--f`` and ``--e``), and employ GHKM algorithm
(``--ghkm``) to perform tree-to-string grammar extraction.
The grammar is constrained by the ``--constrained`` flag, since we
perform exhaustive binarization in the source forest, which may result
in spuriously many minimum synchronous rules. Extracted counts are
stored at `samples/kftt.30k/t2s/model/ghkm-counts`, and models are
stored at `samples/kftt.30k/t2s/model/ghkm-score`. See
`doc/extract.rst` for details.

Features
````````
The extracted grammar at `samples/kftt.30k/t2s/model/ghkm-score`
should be reinterpreted as features as follows:

.. code:: bash

  cd samples/kftt.30k/t2s/model
  ../../../../scripts/cicada-index.py \
	--model-dir . \
	--ghkm \
	--sigtest-inclusive \
	--feature-lexicon \
	--feature-root \
	--feature-internal \
	--feature-height \
	--threads 4

Here, we encode additional features, ``--feature-lexicon`` for lexical
weights, ``--feature-root`` for the generative probabilities,
:math:`\log p(source | root-of-elementary-tree)` and
:math:`\log p(target | root-of-elementary-tree)`.
Also, the number of internal nodes and the height of elementary tree
is used as features (``--feature-internal`` and ``--feature-height``).
See `doc/indexing.rst` for details.

Tuning and testing
``````````````````

During tuning, we use parse forest in hypergraph format as a source
set, which is preprocessed in the similar way as the source side of
the training data. Another difference is the configuration file in
which we use ``--ghkm`` flag:

.. code:: bash

  cd samples/kftt.30k/t2s/tune
  ../../../../scripts/cicada-config.py \
	--tree-grammar ../model/ghkm-index \
	--goal '[x]' \
	--glue '[x]' \
	--fallback-glue \
	--feature-ngram ../../ngram/ngram.5.en.lm \
	--ghkm \
	--beam 1024 > cicada.config

Since the target side is "string", ``[x]`` is the only symbol used in
the target side. Thus, we use ``[x]`` as a goal (``--goal``) and a
glue symbol (``--glue``). The fallback grammar ``--fallback-glue`` use
all the hyperedge in the input hypergraph as a rule, but replaced all
the syntactic labels into ``[x]``.

String-to-tree Model
--------------------



Preprocessing
`````````````

.. code:: bash

  cd samples/kftt.30k/s2t/data
  export stanford=[directory for stanford parser]
  bzcat ../../data/train.en.bz2 | \
  java \
	-mx12g \
	-cp $stanford/stanford-parser.jar:$stanford/stanford-parser-3.2.0-models.jar \
	-tLPP edu.stanford.nlp.parser.lexparser.EnglishTreebankParserParams \
	-tokenized -sentences newline \
	-escaper edu.stanford.nlp.process.PTBEscapingProcessor \
	-outputFormat oneline \
	-outputFormatOptions includePunctuationDependencies \
	edu/stanford/nlp/models/lexparser/englishFactored.ser.gz \
	- | \
  ../../../../progs/cicada_filter_penntreebank \
	--map ../../data/train.en.bz2 \
	--normalize | \
  ../../../../progs/cicada \
	--input-forest \
	--threads 8 \
	--operation binarize:direction=left,order=2 \
	--operation output:no-id=true,file=train.tree.en.gz


Extraction
``````````

.. code:: bash

  cd samples/kftt.30k/s2t/model
  ../../../../scripts/cicada-extract.py \
	--f ../../data/train.ja.bz2 \
	--e ../../data/train.en.bz2 \
	--a ../../alignment/model/aligned.posterior-itg \
	--fe ../data/train.tree.en.gz \
	--model-dir . \
	--lexicon-variational \
	--ghkm \
	--max-scope 2 \
	--threads 4 

Features
````````

.. code:: bash

  cd samples/kftt.30k/s2t/model
  ../../../../scripts/cicada-index.py \
	--model-dir . \
	--ghkm \
	--cky \
	--sigtest-inclusive \
	--feature-lexicon \
	--feature-root \
	--feature-internal \
	--feature-height \
	--threads 4


Tuning and testing
``````````````````

.. code:: bash

  cd samples/kftt.30k/s2t/tune
  ../../../../scripts/cicada-config.py \
	--tree-grammar ../model/ghkm-index \
	--max-span 0 \
	--goal '[ROOT]' \
	--insertion \
	--feature-ngram ../../ngram/ngram.5.en.lm \
	--tree-cky \
	--beam 200 > cicada.config

Tree-to-tree Model
------------------

Preprocessing
`````````````

Extraction
``````````

.. code:: bash

  cd samples/kftt.30k/t2t/tune
  ../../../../scripts/cicada-extract.py \
	--f ../../data/train.ja \
	--e ../../data/train.en \
	--a ../../alignment/model/aligned.posterior-itg \
	--ff ../../t2s/data/train.forest.ja.gz \
	--fe ../../s2t/data/train.tree.en.gz \
	--model-dir . \
	--lexicon-variational \
	--tree \
	--constrained \
	--threads 4

Features
````````

.. code:: bash

  cd samples/kftt.30k/t2t/model
  ../../../../scripts/cicada-index.py \
	--model-dir . \
	--tree \
	--sigtest-inclusive \
	--feature-lexicon \
	--feature-root \
	--feature-internal \
	--feature-height \
	--threads 4

Tuning and testing
``````````````````

.. code:: bash

  cd samples/kftt.30k/t2t/tune
  ../../../../scripts/cicada-config.py \
	--tree-grammar ../model/tree-index \
   	--goal '[ROOT]' \
	--glue '[x]' \
	--fallback-glue \
	--feature-ngram ../../ngram/ngram.5.en.lm \
	--tree \
	--beam 1024 > cicada.config
