===============
cicada-index.py
===============

-----------------------
a tool to index grammar
-----------------------

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-7-25
:Manual section: 1

SYNOPSIS
--------

**cicada-index.py** [*options*]


DESCRIPTIONS
------------

OPTIONS
-------

Input/output
````````````

  --root-dir=DIRECTORY  root directory for outputs
  --model-dir=DIRECTORY
                        model directory (default: ${root_dir}/model)
  --lexical-dir=DIRECTORY
                        lexical transltion table directory (default:
                        ${model_dir})
  --counts-dir=DIRECTORY
                        grammar counts directory (default: ${model_dir})
  --score-dir=DIRECTORY
                        grammar score directory (default: ${model_dir})
  --index-dir=DIRECTORY
                        grammar index directory (default: ${model_dir})
  --temporary-dir=DIRECTORY
                        temporary directory
  --lexicon-source-target=LEXICON
                        lexicon for P(target | source) (default:
                        ${lexical_dir}/lex.f2n)
  --lexicon-target-source=LEXICON
                        lexicon for P(source | target) (default:
                        ${lexical_dir}/lex.n2f)

Models
``````

  --prior=PRIOR         model prior (default: 0.1)
  --feature=FEATURE     feature definitions
  --attribute=ATTRIBUTE
                        attribute definitions
  --phrase              index phrase grammar
  --scfg                index SCFG grammar
  --ghkm                index GHKM (tree-to-string) grammar
  --tree                index tree-to-tree grammar
  --cky                 CKY mode indexing for tree-grammar
  --reordering          reordering for phrase grammar

Additional features
```````````````````

  --feature-root        generative probability
  --feature-fisher      Fisher's exact test
  --feature-type        observation probability
  --feature-singleton   singleton features
  --feature-cross       cross features
  --feature-unaligned   unaligned features
  --feature-internal    internal features
  --feature-height      height features
  --feature-lexicon     Lexical features
  --feature-model1      Model1 features
  --feature-noisy-or    noisy-or features
  --feature-insertion-deletion
                        insertion/deletion features
  --prefix-feature=PREFIX_FEATURE
                        feature name prefix (default: )
  --prefix-attribute=PREFIX_ATTRIBUTE
                        attribute name prefix (default: )
  --threshold-insertion=THRESHOLD_INSERTION
                        threshold for insertion (default: 0.5)
  --threshold-deletion=THRESHOLD_DELETION
                        threshold for deletion (default: 0.5)
  --quantize            perform quantization


Pruning
```````

  --kbest=KBEST         kbest max count of rules (default: 0)
  --kbest-count         kbest option: use count for sorting
  --kbest-joint         kbest option: use joint probability for sorting
  --kbest-source        kbest option: use source probability (P(f|e)) for
                        sorting
  --kbest-target        kbest option: use target probability (P(e|f)) for
                        sorting
  --cutoff=CUTOFF       cutoff count of rules (default: 0)
  --threshold=THRESHOLD
                        probability threshold of rules (default: 0)
  --sigtest=SIGTEST     significance testing threshold relative to 1-1-1-N
                        log-p-value (or \epsilon in "discarding most of the
                        phrasetable") (default: 0)
  --sigtest-inclusive   significance testing which includes 1-1-1-N event
                        (this will assign --sigtest -0.001)
  --sigtest-exclusive   significance testing which excludes 1-1-1-N event
                        (this will assign --sigtest +0.001)

Others
``````

  --max-malloc=MALLOC   maximum memory in GB (default: 8)
  --cicada-dir=DIRECTORY
                        cicada directory
  --mpi-dir=DIRECTORY   MPI directory
  --threads=THREADS     # of thrads for thread-based parallel processing
  --mpi=MPI             # of processes for MPI-based parallel processing.
                        Identical to --np for mpirun
  --mpi-host=HOSTS      list of hosts to run job. Identical to --host for
                        mpirun
  --mpi-host-file=FILE  host list file to run job. Identical to --hostfile for
                        mpirun
  --mpi-options=OPTION  additional MPI options
  --pbs                 PBS for launching processes
  --pbs-queue=NAME      PBS queue for launching processes (default: ltg)
  --debug=DEBUG         debug level
  -h, --help            show this help message and exit

EXAMPLES
--------


SEE ALSO
--------
