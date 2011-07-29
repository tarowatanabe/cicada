#!/bin/sh

me=`basename $0`

root=""
cicada=/data/lttools/decoder/cicada/bin
openmpi=""

tstset=""
refset=""

## # of processes, # of cores
np=1
nc=2
hosts=""
hosts_file=""

### decoding config
config=""

### linear learning
iteration=10
weights_init=""
C=1e-3

### qsubs
mem_single=1gb
mem_mpi=8gb
queue=ltg

exit_missing_arg="\
echo \"$me: option \\\`\$1' requires an argument\" >&2
echo \"\$help\" >&2
exit 1"

usage="\
$me [options]
  General options
  --root                    root directory
  --cicada                  cicada directory (required)
  --mpi                     MPI directory
  --host, --hosts           MPI hosts
  --hostfile, --host-file   MPI host file
  -q, --queue               PBS queue
  -n, --num                 # of processes to run
  
  Decoding options
  -c, --config              Configuration file (required)
  
  Training options
  -i, --iteration           MERT iterations
  -w, --weights             initial weights
  -C, --C                   hyperparameter
  -t, --test, --tstset      tuning data (required)
  -r, --reference, --refset reference translations (required)

  -h, --help                help message
"

while test $# -gt 0 ; do
  case $1 in
  --root )
    test $# = 1 && eval "$exit_missing_arg"
    root=$2
    shift; shift ;;
  --cicada )
    test $# = 1 && eval "$exit_missing_arg"
    cicada=$2
    shift; shift ;;
  --mpi )
    test $# = 1 && eval "$exit_missing_arg"
    openmpi=$2
    shift; shift ;;
  --host | --hosts )
    test $# = 1 && eval "$exit_missing_arg"
    hosts=$2
    shift; shift ;;
  --hostfile | --host-file )
    test $# = 1 && eval "$exit_missing_arg"
    host_file=$2
    shift; shift ;;
  --queue | -q )
    test $# = 1 && eval "$exit_missing_arg"
    queue=$2
    shift; shift ;;
  --num | -n )
    test $# = 1 && eval "$exit_missing_arg"
    np=$2
    shift; shift ;;

  ## training
  --iteration | -i )
    test $# = 1 && eval "$exit_missing_arg"
    iteration=$2
    shift; shift ;;
  --weights | -w )
    test $# = 1 && eval "$exit_missing_arg"
    weights_init=$2
    shift; shift ;;
  --C | -C )
    test $# = 1 && eval "$exit_missing_arg"
    C=$2
    shift; shift ;;

  --config | -c )
    test $# = 1 && eval "$exit_missing_arg"
    config=$2
    shift; shift ;;

### test set and reference set
  --test | -t | --tstset )
    test $# = 1 && eval "$exit_missing_arg"
    tstset=$2
    shift; shift ;;
  --reference | -r | --refset )
    test $# = 1 && eval "$exit_missing_arg"
    refset=$2
    shift; shift ;;
  --help | -h )
    echo "$usage"
    exit ;;
### error...
   -* )
    exec >&2
    echo "$me: invalid option $1"
    echo "$help"
    exit 1 ;;
  * )
    break ;;
  esac
done

if test "$tstset" = ""; then
  echo "specify development data"
  exit -1
fi
if test "$refset" = ""; then
  echo "specify reference data"
  exit -1
fi
if test "$config" = ""; then
  echo "specify config file"
  exit -1
fi

if test "$weights_init" != ""; then
  if test ! -e $weights_init; then
    echo "no initial weights: $weights_init ?"
    exit -1
  fi
fi

if test "$openmpi" != ""; then
  openmpi="${openmpi}/"
fi

if test "$root" != ""; then
  root="${root}/"
fi

### working dir..
workingdir=`pwd`

### this is a test, whether we can run under cluster or not...
qsub=`which qsub 2> /dev/null`

mpinp=""
if test "$qsub" = ""; then
  mpinp="--np $np"
  if test "$hosts" != ""; then
    mpinp="$mpinp --host $hosts"
  fi
  if test "$host_file" != ""; then
    mipnp="$mpinp --hostfile $host_file"
  fi
fi


qsubsingle() {
  name=$1
  shift

  logfile=/dev/null
  if [ "$1" = "-l" ]; then
    logfile=$workingdir/$2
    shift 2
  fi

  if test "$qsub" != ""; then
    (
      echo "#!/bin/sh"
      echo "#PBS -N $name"
      echo "#PBS -W block=true"
      echo "#PBS -e $logfile"
      echo "#PBS -o /dev/null"
      echo "#PBS -q $queue"
      echo "#PBS -l select=1:ncpus=1:mem=${mem_singl}"
      echo "cd $workingdir"
      echo "$@"
    ) |
    qsub -S /bin/sh > /dev/null
  else
    $@ >& $logfile
  fi
}

qsubwrapper() {
  name=$1
  shift

  logfile=/dev/null
  if [ "$1" = "-l" ]; then
    logfile=$workingdir/$2
    shift 2
  fi

  if test "$qsub" != ""; then
    (
      echo "#!/bin/sh"
      echo "#PBS -N $name"
      echo "#PBS -W block=true"
      echo "#PBS -e $logfile"
      echo "#PBS -o /dev/null"
      echo "#PBS -q $queue"
      echo "#PBS -l select=$np:ncpus=$nc:mpiprocs=$nc:mem=${mem_mpi}"
      echo "#PBS -l place=scatter"
      echo "cd $workingdir"
      echo "$@"
    ) |
    qsub -S /bin/sh > /dev/null
  else
    $@ >& $logfile
  fi
}

for ((iter=1;iter<=iteration; ++ iter)); do
  echo "iteration: $iter"
  iter_prev=`expr $iter - 1`

  #
  # as our initial weights, use the ${weights_init} if exists
  #
  weights="weights-one=true" 
  if test "${weights_init}" != ""; then
    weights="weights=${weights_init}"  
  fi

  if test -e "${root}weights.$iter_prev"; then
    weights="weights=${root}weights.$iter_prev"
  fi

  # generate translation and generate hypergraph (kbest=0)
  # instead of compose-cky, we use parse-cky with ${weights_init} learned by learn.sh!

  qsubsingle config ${cicada}/cicada_filter_config \
      --weights $weights \
      --directory ${root}kbest-$iter \
      --input $config \
      --output ${root}cicada.config.$iter
 
#  qsubwrapper decode -l ${root}decode.$iter.log ${openmpi}mpirun $mpinp $cicada/cicada_mpi \
   qsubsingle decode -l ${root}decode.$iter.log $cicada/cicada \
	--input $tstset \
	--config ${root}cicada.config.$iter \
	\
	--debug


  qsubsingle eval $cicada/cicada_eval --refset $refset --tstset ${root}kbest-$iter --output ${root}eval-$iter.1best --scorer bleu:order=4

  ### compute oracles
#  qsubwrapper oracle -l oracle.$iter.log $openmpi/mpirun $mpinp $cicada/cicada_oracle_kbest_mpi \
  qsubsingle oracle -l ${root}oracle.$iter.log $cicada/cicada_oracle_kbest \
        --refset $refset \
        --tstset ${root}kbest-$iter \
        --output ${root}kbest-${iter}.oracle \
        --directory \
        --scorer  bleu:order=4,exact=true \
        \
        --debug


  ### kbests upto now...
  tsts=""
  orcs=""
  for ((i=1;i<=$iter;++i)); do
    if test -e ${root}kbest-$i; then
      tsts="$tsts ${root}kbest-$i"
      orcs="$orcs ${root}kbest-${i}.oracle"
    fi
  done

  qsubsingle learn -l ${root}learn.$iter.log $cicada/cicada_learn_kbest \
                        --kbest  $tsts \
                        --oracle $orcs \
                        --output ${root}weights.$iter \
                        \
                        --learn-linear \
                        --solver 1 \
                        --C $C \
                        \
                        --debug=2


done 
