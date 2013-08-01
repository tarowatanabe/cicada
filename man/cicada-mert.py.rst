==============
cicada-mert.py
==============

----------------------------------
a tool to learn parameters by MERT
----------------------------------

:Author: Taro Watanabe <taro.watanabe@nict.go.jp>
:Date:   2013-8-1
:Manual section: 1

SYNOPSIS
--------

**cicada-mert.py** [*options*]

OPTIONS
-------

  --root-dir=DIRECTORY  root directory for outputs
  --prefix=PREFIX       prefix for outputs (default: mert)
  --srcset=FILE         training data
  --refset=FILE         reference translations
  --config=CONFIG       cicada config file
  --iteration=ITERATION
                        # of iterations (default: 10)
  --iteration-first=ITERATION
                        The first iteration (default: 1)
  --weights=FILE        initial weights
  --bound-lower=FILE    lower bounds for weights
  --bound-upper=FILE    upper bounds for weights
  --parameter-lower=PARAMETER_LOWER
                        lower parameter value (default: -1.0)
  --parameter-upper=PARAMETER_UPPER
                        upper parameter value (default: 1.0)
  --mert-options=MERT_OPTIONS
                        other MERT options
  --direction=DIRECTION
                        # of random directions (default: 8)
  --restart=RESTART     # of random restarts (default: 2)
  --scorer=SCORER       scorer for oracle computation (default:
                        bleu:order=4,exact=true)
  --kbest=KBEST         kbest size (default: 0)
  --forest              forest based learning
  --iterative           perform iterative learning
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
