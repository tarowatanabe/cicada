Indexer for Phrase/Rule/Tree-fragment


Usage: cicada-index.py [options]

Options:
  --root-dir=DIRECTORY  root directory for outputs
  --model-dir=DIRECTORY
                        model directory (default: ${root_dir}/model)
  --lexical-dir=DIRECTORY
                        lexical transltion table directory (default:
                        ${model_dir})
  --lexical-source-target=LEXICON
                        lexicon for P(target | source) (default:
                        ${lexical_dir}/lex.f2n)
  --lexical-target-source=LEXICON
                        lexicon for P(source | target) (default:
                        ${lexical_dir}/lex.n2f)
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
  --feature-root        generative probability
  --feature-type        observation probability
  --feature-singleton   singleton features
  --feature-cross       cross features
  --feature-unaligned   unaligned features
  --feature-lexicon     compute Model1 features
  --feature-model1      compute Model1 features
  --feature-noisy-or    compute noisy-or features
  --feature-insertion-deletion
                        compute insertion/deletion features
  --threshold-insertion=THRESHOLD_INSERTION
                        threshold for insertion (default: 0.01)
  --threshold-deletion=THRESHOLD_DELETION
                        threshold for deletion (default: 0.01)
  --quantize            perform quantization
  --kbest=KBEST         kbest max count of rules
  --max-malloc=MALLOC   maximum memory in GB (default: 4)
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
  --pbs                 PBS for launching processes
  --pbs-queue=NAME      PBS queue for launching processes (default: ltg)
  --debug=DEBUG         
  -h, --help            show this help message and exit

