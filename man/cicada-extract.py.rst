=================
cicada-extract.py
=================

---------------------------------------------
a tool to extract grammar from bilingual data
---------------------------------------------

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-2-28
:Manual section: 1

SYNOPSIS
--------

**cicada-extract.py** [*options*]

OPTIONS
-------

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
                        alignment methods (default: grow-diag-final-and)
  --first-step=STEP     first step (default: 4)
  --last-step=STEP      last step  (default: 6)
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

EXAMPLES
--------




SEE ALSO
--------
