===================
cicada-alignment.py
===================

-------------------------------------------------
a tool to learn word alignment for bilingual data
-------------------------------------------------

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-8-1
:Manual section: 1

SYNOPSIS
--------

**cicada-alignment.py** [*options*]

DESCRIPTION
-----------

Perform word alignment given a bilingual data. Word alignment training
is performed in three steps:

1. Source and target words are clustered, and clustering
   results are placed in the `corpus` directory (or specified by
   ``--corpus-dir`` option).
2. Word alignment models are learned and Viterbi alignment is
   computed, then, all the parameters are put on two directories,
   `giza.src-trg`  and `giza.trg-src` directories (or specified by
   ``--giza-f2e`` and ``--giza-e2f`` options).
3. Two models in opposite directions are combined, and put in `model`
   directory (or specified by ``--model-dir``).

Internally, it calls `cicada_alignment_model1(1)`,
`cicada_alignment_hmm(1)` or `cicada_alignment_model4(1)` for word
alignment learning and Viterbi alignment computation depending on the
final word alignment model. In addition, `cicada_alignment(1)` may be
used when combining word alignment using a heuristic,
i.e. "grow-diag-final-and".

OPTIONS
-------

Input/output
````````````

The bilingual data is provided in two ways. A standard way is to
directly specify the files using the options ``--f`` and ``--e``. An
alternative is to use the ``--corpus`` option which is a prefix of
training data, and ``--f`` and ``--e`` are "suffix" for the
prefix. For instance, if your training data is located at ``data``
directory with shared prefix:
::

  data/prefix.f
  data/prefix.e

then, you can specify ``--corpus data/prefix --f f --e e``, or
``--f data/prefix.f --e data/prefix.e``.

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
                        Note that this is not the final word alignment
			produced by the **cicada-alignment.py**, but
			initial word alignment which are used as
			constraints for alignment model training.
  --sf=FILE-OR-SUFFIX   source (or 'French')  span file or suffix
  --se=FILE-OR-SUFFIX   target (or 'English') span file or suffix
  --ff=FILE-OR-SUFFIX   source (or 'French')  forest file or suffix
  --fe=FILE-OR-SUFFIX   target (or 'English') forest file or suffix
  --corpus=CORPUS       bilingual trainging corpus prefix

Learning options
````````````````

  --alignment=ALIGNMENT
                        alignment methods (default: posterior-itg)
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

Smoothing
`````````

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

Others
``````

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

`cicada_alignment(1)`,
`cicada_alignment_model1(1)`,
`cicada_alignment_hmm(1)`,
`cicada_alignment_model4(1)`
