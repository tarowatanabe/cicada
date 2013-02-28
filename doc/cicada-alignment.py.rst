===================
cicada-alignment.py
===================

-------------------------------------------------
a tool to learn word alignment for bilingual data
-------------------------------------------------

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-2-28
:Manual section: 1

SYNOPSIS
--------

**cicada-alignment.py** [*options*]

OPTIONS
-------

  --root-dir=DIRECTORY  root directory for outputs
  --corpus-dir=PREFIX   corpus directory (default: ${root_dir}/corpus)
  --giza-f2e=DIRECTORY  giza directory for P(f|e) (default:
                        ${root_dir}/giza.${f}-${e})
  --giza-e2f=DIRECTORY  giza directory for P(e|f) (default:
                        ${root_dir}/giza.${e}-${f})
  --model-dir=DIRECTORY
                        model directory (default: ${root_dir}/model)
  --alignment-dir=DIRECTORY
                        alignment directory (default: ${model_dir})
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
  --first-step=STEP     first step (default: 1)
  --last-step=STEP      last step  (default: 3)
  --iteration-cluster=ITERATION
                        word cluter iterations (default: 50)
  --iteration-model1=ITERATION
                        Model1 iteratins (default: 5)
  --iteration-hmm=ITERATION
                        HMM iteratins    (default: 5)
  --iteration-model4=ITERATION
                        Model4 iteratins    (default: 5)
  --cluster=CLUSTER     # of clusters (default: 50)
  --symmetric           symmetric training
  --posterior           posterior constrained training
  --dynamic             dynamically recompute base alignment
  --variational         variational Bayes estimates
  --l0                  L0 regularization
  --p0=P0               parameter for NULL alignment (default: 0.01)
  --insertion-p1=P1     parameter for NULL insertion (default: 0.01)
  --prior-lexicon=PRIOR
                        lexicon model prior (default: 0.01)
  --prior-alignment=PRIOR
                        alignment model prior (default: 0.01)
  --prior-distortion=PRIOR
                        distortion model prior (default: 0.01)
  --prior-fertility=PRIOR
                        fertility model prior (default: 0.01)
  --smooth-lexicon=SMOOTH
                        lower-bound parameter for lexicon model (default:
                        1e-100)
  --smooth-alignment=SMOOTH
                        lower-bound parameter for alignment model (default:
                        1e-100)
  --smooth-distortion=SMOOTH
                        lower-bound parameter for distortion model (default:
                        1e-100)
  --smooth-fertility=SMOOTH
                        lower-bound parameter for fertility model (default:
                        1e-100)
  --l0-alpha=L0_ALPHA   L0 regularization parameter (default: 100)
  --l0-beta=L0_BETA     L0 regularization parameter (default: 0.01)
  --cicada-dir=DIRECTORY
                        cicada directory
  --threads=THREADS     # of thrads for thread-based parallel processing
  --max-malloc=MALLOC   maximum memory in GB (default: 8)
  --pbs                 PBS for launching processes
  --pbs-queue=NAME      PBS queue for launching processes (default: ltg)
  --debug=DEBUG         debug level
  -h, --help            show this help message and exit


EXAMPLES
--------




SEE ALSO
--------
