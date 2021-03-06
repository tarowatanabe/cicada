===============
cicada-learn.py
===============

--------------------------
a tool to learn parameters
--------------------------

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-8-1
:Manual section: 1

SYNOPSIS
--------

**cicada-learn.py** [*options*]

OPTIONS
-------

  --root-dir=DIRECTORY  root directory for outputs
  --prefix=PREFIX       prefix for outputs (default: learn)
  --srcset=FILE         training data
  --refset=FILE         reference translations
  --config=CONFIG       cicada config file
  --iteration=ITERATION
                        # of iterations (default: 10)
  --iteration-first=ITERATION
                        The first iteration (default: 1)
  --weights=FILE        initial weights
  --regularize-l1=REGULARIZER
                        L1 regularization
  --regularize-l2=REGULARIZER
                        L2 regularization
  --regularize-lambda=REGULARIZER
                        regularization hyperparameter
  --regularize-oscar=REGULARIZER
                        OSCAR regularization
  --scorer=SCORER       scorer for oracle computation (default:
                        bleu:order=4,exact=true)
  --scorer-cube=SIZE    cube size for oracle computation (default: 400)
  --learn=LEARN         learning algorithms from [softmax, svm, linear,
                        pegasos, mira, cw, arow, nherd, cp, mcp, xbleu]
                        (default: xbleu)
  --learn-options=OPTION
                        additional learning options
  --kbest=KBEST         kbest size (default: 0)
  --forest              forest based learning
  --merge               perform kbest merging
  --interpolate=INTERPOLATE
                        perform weights interpolation
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
  --mpi-options=OPTION  additional MPI options
  --threads=THREADS     # of thrads for thread-based parallel processing
  --pbs                 PBS for launching processes
  --pbs-queue=NAME      PBS queue for launching processes (default: ltg)
  --debug=DEBUG         debug level
  -h, --help            show this help message and exit


EXAMPLES
--------


SEE ALSO
--------
