=================
cicada-extract.py
=================

---------------------------------------------
a tool to extract grammar from bilingual data
---------------------------------------------

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-7-31
:Manual section: 1

SYNOPSIS
--------

**cicada-extract.py** [*options*]


DESCRIPTIONS
------------

OPTIONS
-------

Input/output
````````````

  --root-dir=DIRECTORY  root directory for outputs
  --corpus-dir=PREFIX   corpus directory (default: ${root_dir}/corpus)
  --model-dir=DIRECTORY
                        model directory (default: ${root_dir}/model)
  --alignment-dir=DIRECTORY
                        alignment directory (default: ${model_dir})
  --lexical-dir=DIRECTORY
                        lexical transltion table directory (default:
                        ${model_dir})
  --counts-dir=DIRECTORY
                        grammar counts directory (default: ${model_dir})
  --score-dir=DIRECTORY
                        grammar score directory (default: ${model_dir})
  --temporary-dir=DIRECTORY
                        temporary directory
  --f=FILE-OR-SUFFIX    source (or 'French')  language file or suffix
  --e=FILE-OR-SUFFIX    target (or 'English') language file or suffix
  --a=FILE-OR-SUFFIX    alignment file or suffix
  --sf=FILE-OR-SUFFIX   source (or 'French')  span file or suffix
  --se=FILE-OR-SUFFIX   target (or 'English') span file or suffix
  --ff=FILE-OR-SUFFIX   source (or 'French')  forest file or suffix
  --fe=FILE-OR-SUFFIX   target (or 'English') forest file or suffix
  --corpus=CORPUS       bilingual trainging corpus prefix
  --alignment=ALIGNMENT
                        alignment methods (default: posterior-itg)

Extraction steps
````````````````

  --first-step=STEP     first step (default: 1)
  --last-step=STEP      last step  (default: 3)

Lexicon
```````

  --lexicon-inverse     use inverse alignment
  --lexicon-prior=PRIOR
                        lexicon model prior (default: 0.1)
  --lexicon-variational
                        variational Bayes estimates
  --lexicon-l0          L0 regularization
  --lexicon-l0-alpha=LEXICON_L0_ALPHA
                        L0 regularization parameter (default: 100)
  --lexicon-l0-beta=LEXICON_L0_BETA
                        L0 regularization parameter (default: 0.01)

Models
``````

  --phrase              extract phrase
  --scfg                extract SCFG
  --ghkm                extract GHKM (tree-to-string)
  --tree                extract tree-to-tree
  --non-terminal=NON_TERMINAL
                        default non-terminal for GHKM rule (default: [x])


  --max-sentence-length=LENGTH
                        maximum sentence size (default: 0 == no limit)
  --max-span-source=LENGTH
                        maximum source span size (default: 15)
  --max-span-target=LENGTH
                        maximum target span size (default: 15)
  --min-hole-source=LENGTH
                        minimum source hole size (default: 1)
  --min-hole-target=LENGTH
                        minimum target hole size (default: 1)
  --max-length=LENGTH   maximum terminal length (default: 7)
  --max-fertility=FERTILITY
                        maximum terminal fertility (default: 4)
  --max-nodes=NODES     maximum internal nodes (default: 7)
  --max-height=HEIGHT   maximum rule height (default: 4)
  --max-compose=COMPOSE
                        maximum rule composition (default: 0)
  --max-rank=RANK       maximum rule rank (default: 2)
  --max-scope=SCOPE     maximum rule scope (default: 0)
  --cutoff=CUTOFF       cutoff counts (default: 0.0)
  --collapse-source     collapse source side for CKY parsing
  --collapse-target     collapse target side for CKY parsing
  --exhaustive          exhaustive extraction in SCFG, GHKM and Tree
  --constrained         constrained extraction in SCFG, GHKM and Tree
  --project             project non-terminal symbols in GHKM
  --sentential          extract sentential rule

Others
``````

  --max-malloc=MALLOC   maximum memory in GB (default: 8)
  --cicada-dir=DIRECTORY
                        cicada directory
  --mpi-dir=DIRECTORY   MPI directory
  --mpi=MPI             # of processes for MPI-based parallel processing.
                        Identical to --np for mpirun
  --mpi-host=HOSTS      list of hosts to run job. Identical to --host for
                        mpirun
  --mpi-host-file=FILE  host list file to run job. Identical to --hostfile for
                        mpirun
  --threads=THREADS     # of thrads for thread-based parallel processing
  --pbs                 PBS for launching processes
  --pbs-queue=NAME      PBS queue for launching processes (default: ltg)
  --debug=DEBUG         debug level
  -h, --help            show this help message and exit

DETAILS
-------

Internally, cicada-extract.py calls following binaries:

cicada_extract_phrase{,_mpi}
	Extract phrase
	Output scores:

	       lhs ||| rhs ||| count(lhs, rhs) \
	       	   count(prev, mono, lhs, rhs) \
		   count(prev, swap, lhs, rhs) \
		   count(next, mono, lhs, rhs) \
		   count(next, swap, lhs, rhs)

cicada_extract_rule{,_mpi}
	Extract synchronous-CFG + syntax augmentation (aka SAMT) when extracted with "span" data.
	You can generate span by "cicada_filter_penntreebank"
	Output scores:

	       root lhs ||| root rhs ||| count(lhs, rhs)

cicada_extract_ghkm{,_mpi}
	Extract tree-to-string or string-to-tree rules by GHKM
	Output scores:

	       lhs-xRS ||| rhs-xRS ||| count(lhs, rhs)

cicada_extract_tree{,_mpi}
	Extract tree-to-tree rules by GHKM
	Output scores:

::

	       lhs-xRS ||| rhs-xRS ||| count(lhs, rhs)

After counts collection, you can summarize them by cicada_extract_counts{,_mpi}
Which will output:

lhs ||| rhs ||| alignments ||| counts(lhs, rhs) ||| counts(lhs) ||| counts(rhs) ||| observed(lhs) observed(rhs) 

In addtion, we will dump root-joint.gz, root-source.gz and root-target.gz. root-joint.gz contains:

root(lhs)root(rhs) ||| counts(root(lhs)root(rsh)) ||| observed(lhs, rhs)

while root-source.gz looks like:

root(lhs) ||| counts(root(lhs)) ||| observed(lhs)

while root-target.gz looks like:

root(rhs) ||| counts(root(rhs)) ||| observed(rhs)

You can easily transform the counts into probabilities by maximum likelihood estimates, or use observed counts
to perform Dirichlet prior smoothing (default) by running "cicada-index.py".
The cicada-index.py transforms collected counts into a set of feature values, then, encodes the grammar into a binary format.
Internally, the indexer calls:

   cicada_filter_extract:
	extract only nbet of target variation for each source side,
   	measured by its joint frequency of lhs and rhs.

   cicada_filter_extract_phrase:
	Dump in moses or cicada format. Also, you can dump
	lexicalied reordering table.
   				 
   cicada_filter_extract_scfg:
	Dump in cicada format for synchronous-CFG. You can also 
	add features for lhs given root(lhs) and rhs given root(rhs)

   cicada_filter_extract_ghkm:
	Dump in cicada format for tree-to-string, string-to-tree, tree-to-tree.

and calls
    cicada_index_grammar
	for indexing phrase/scfg
    
    or

    cicada_index_tree_grammar
	for indexing tree-to-string, string-to-tree, tree-to-tree rules


EXAMPLES
--------




SEE ALSO
--------
