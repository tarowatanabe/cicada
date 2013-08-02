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
:math:`\log p(source | root(source))` and
:math:`\log p(target | root(target))`.
Also, the number of internal nodes and the height of elementary tree
is used as features (``--feature-internal`` and ``--feature-height``).
See `doc/indexing.rst` for details.

Tuning and testing
``````````````````

During tuning, we use parse forest in hypergraph format as a source
set, which is preprocessed in the similar way as the source side of
the training data. Another difference is the configuration file in
which we use ``--tree`` flag:

.. code:: bash

  cd samples/kftt.30k/t2s/tune
  ../../../../scripts/cicada-config.py \
	--tree-grammar ../model/ghkm-index \
	--goal '[x]' \
	--glue '[x]' \
	--fallback-glue \
	--feature-ngram ../../ngram/ngram.5.en.lm \
	--tree \
	--beam 1024 > cicada.config

Since the target side is "string", ``[x]`` is the only symbol used in
the target side. Thus, we use ``[x]`` as a goal (``--goal``) and a
glue symbol (``--glue``). The fallback grammar ``--fallback-glue`` use
all the hyperedge in the input hypergraph as a rule, but replaced all
the syntactic labels into ``[x]``.

String-to-tree Model
--------------------

In a string-to-tree model, input is a string, or a sentence, which is
parsed by the source side of synchronous-TSG. Translation forest is
constructed by the projected elementary trees.  Note that
the string-to-tree model is very slow in practice since we need to
parse inputs using the source side of the synchronous grammar. In
addition, it may consume large memory especially when the synchronous
grammar is ambiguous.

Preprocessing
`````````````

In this example, we use the Stanford Parser (http://nlp.stanford.edu/software/lex-parser.shtml), 
for parsing the English side of bilingual data, but any parser will work.

.. code:: bash

  cd samples/kftt.30k/s2t/data
  SP=[directory for stanford parser]
  bzcat ../../data/train.en.bz2 | \
  java \
	-mx12g \
	-cp $SP/stanford-parser.jar:$SP/stanford-parser-3.2.0-models.jar \
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

Here, we assume that English sentences are tokenized, but further
escaped to match with the Penntreebank standard, like ``(`` into
``-LRB-`` etc. (see ``edu.stanford.nlp.process.PTBEscapingProcessor``).
Thus, the penntreebank to hypergraph converter, ``cicada_filter_penntreebank``
re-map the escaped terminal symbols again via ``--map`` argument.
The constituency labels are also normalized (``--normalize``) so that
we can use ``COMMA`` as a label for ``,``. The hypergraph is further
binarized in a left-heavy direction (``binarize:direction=left,order=2``)
and preserves only 2 non-terminal history for the binarized symbols.

Extraction
``````````

We use the word alignment extracted in `samples/kftt.30k/alignment`
and extract synchronous rules:

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

In stead of using the source side forest in the above tree-to-string
model, we use a tree in the target side (``--fe``). Since the tree
structure is not exhaustively binarized, we do not constrained as in
the previous example, but constrained the extracted grammar so that
the maximum scope is 2 (``--max-scope 2``) which greatly affect the
parsing speed.

Features
````````

The extracted grammar at `samples/kftt.30k/st2/model/ghkm-score`
should be reinterpreted as features as follows:

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

Note that we have ``--cky`` flag which indicates that the indexed
model is suitable for parsing a string by the CKY algorithm.

Tuning and testing
``````````````````

For tuning and testing, input is a sentence, as in SCFG. We use
``--tree-cky`` as an algorithm to parse inputs with the learned
synchronous grammar. In addition, we add insertion grammar as in SCFG
which can handle OOVs.

.. code:: bash

  cd samples/kftt.30k/s2t/tune
  ../../../../scripts/cicada-config.py \
	--tree-grammar ../model/ghkm-index \
	--max-span 10 \
	--goal '[ROOT]' \
	--tree-straight \
	--insertion \
	--feature-ngram ../../ngram/ngram.5.en.lm \
	--tree-cky \
	--beam 100 > cicada.config

In this example, we add straight glue rules as in SCFG, but using the
STSG (``--tree-straight``), but limit the maximum span to 10 and beam
size to 100 for practical reason.

Tree-to-tree Model
------------------

An easy extension to the tree-to-string model is to incorporate target
side syntax into the synchronous grammar. In this example, we use the
Japanese trees from the tree-to-string example and the English trees
from the string-to-tree example and create a model.

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
	--exhaustive \
	--threads 4

Here, we use the ``--tree`` flag to indicate the tree-to-tree model
extraction using the source and target syntactic trees (``--ff`` and
``--fe``). Also, we extract rules exhaustively (``--exhaustive`` flag)
by considering all possible attachments of non aligned words.

Features
````````

The extracted model is interpreted using the tree-to-tree model (``--tree``).

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

Similar to the tree-to-string model we use the same tree composition
algorithm (``--tree``) but ``[ROOT]`` as a goal symbol which is the
syntactic label in the target trees. As a glue rules, we use ``[x]``
as described in the tree-to-string model.

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
