#!/bin/sh

me=`basename $0`

root=""
cicada=""
openmpi=""

devset=""
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
merge="no"
merge_ratio=0.0

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
  -d, --dev, --devset       tuning data (required)
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
  --merge )
    merge=yes
    shift ;;
  --merge-ratio )
    test $# = 1 && eval "$exit_missing_arg"
    merge_ratio=$2
    shift; shift ;;

  --config | -c )
    test $# = 1 && eval "$exit_missing_arg"
    config=$2
    shift; shift ;;

### test set and reference set
  --dev | -d | --devset )
    test $# = 1 && eval "$exit_missing_arg"
    devset=$2
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

if test "$devset" = ""; then
  echo "specify development data"
  exit 1
fi
if test "$refset" = ""; then
  echo "specify reference data"
  exit 1
fi
if test "$config" = ""; then
  echo "specify config file"
  exit 1
fi
if test "$cicada" = ""; then
  echo "no cicada dir?"
  exit 1
fi

if test "$weights_init" != ""; then
  if test ! -e $weights_init; then
    echo "no initial weights: $weights_init ?"
    exit 1
  fi
fi

if test "$openmpi" != ""; then
  openmpi=`echo "${openmpi}/" | sed -e 's/\/\/$/\//'`
fi

if test "$root" != ""; then
 root=`echo "${root}/" | sed -e 's/\/\/$/\//'`
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

qsubwrapper() {
  name=$1
  shift
  
  logfile=""
  while test $# -gt 0 ; do
  case $1 in
  -l )
    test $# = 1 && eval "$exit_missing_arg"
    logfile=$2
    shift; shift ;;
  -* )
    exec >&2
    echo "$me: invalid option $1"
    exit 1 ;;
  * )
    break ;;
  esac
  done

  stripped=`expr "$1" : '^\(.*\)_mpi$'`
  if test "$stripped" = ""; then
    stripped=$1
  fi

  if test "$qsub" != ""; then
    (
      echo "#!/bin/sh"
      echo "#PBS -N $name"
      echo "#PBS -W block=true"
      echo "#PBS -e /dev/null"
      echo "#PBS -o /dev/null"
      echo "#PBS -q $queue"
      if test "$stripped" != "$1" -a $np -gt 1; then
        echo "#PBS -l select=$np:ncpus=$nc:mpiprocs=$nc:mem=${mem_mpi}"
        echo "#PBS -l place=scatter"
      else
        echo "#PBS -l select=1:ncpus=1:mem=${mem_single}"
      fi
      echo "cd $workingdir"

      if test "$stripped" != "$1" -a $np -gt 1; then
        if test "$logfile" != ""; then
          echo "${openmpi}mpirun $mpinp $@ >& $logfile"
        else
          echo "${openmpi}mpirun $mpinp $@"
        fi
      else
	## shift here!
	shift;
	if test "$logfile" != ""; then
          echo "$stripped $@ >& $logfile"
        else
          echo "$stripped $@"
        fi
      fi
    ) |
    qsub -S /bin/sh || exit 1
  else
    if test "$stripped" != "$1" -a $np -gt 1; then
      if test "$logfile" != ""; then
        ${openmpi}mpirun $mpinp $@ >& $logfile || exit 1
      else
        ${openmpi}mpirun $mpinp $@ || exit 1
      fi
    else
      shift
      if test "$logfile" != ""; then
        $stripped $@ >& $logfile || exit 1
      else
        $stripped $@ || exit 1
      fi
    fi
  fi
}

for ((iter=1;iter<=iteration; ++ iter)); do
  echo "iteration: $iter" >&2
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

  ### setup config file
  echo "generate config file ${root}cicada.config.$iter" >&2
  qsubwrapper config ${cicada}/cicada_filter_config \
      --weights $weights \
      --directory ${root}kbest-$iter \
      --input $config \
      --output ${root}cicada.config.$iter || exit 1
  
  ### actual decoding
  echo "decoding ${root}kbest-$iter" >&2
  qsubwrapper decode -l ${root}decode.$iter.log $cicada/cicada_mpi \
	--input $devset \
	--config ${root}cicada.config.$iter \
	\
	--debug || exit 1

  ### BLEU
  echo "BLEU ${root}eval-$iter.1best" >&2
  qsubwrapper eval $cicada/cicada_eval \
      --refset $refset \
      --tstset ${root}kbest-$iter \
      --output ${root}eval-$iter.1best \
      --scorer bleu:order=4 || exit 1

  ### compute oracles
  echo "oracle translations ${root}kbest-${iter}.oracle" >&2
  qsubwrapper oracle -l ${root}oracle.$iter.log $cicada/cicada_oracle_kbest_mpi \
        --refset $refset \
        --tstset ${root}kbest-$iter \
        --output ${root}kbest-${iter}.oracle \
        --directory \
        --scorer  bleu:order=4,exact=true \
        \
        --debug || exit 1


  ### kbests upto now...
  tstset=""
  orcset=""
  for ((i=1;i<=$iter;++i)); do
    if test -e ${root}kbest-$i; then
      tstset="$tstset ${root}kbest-$i"
      orcset="$orcset ${root}kbest-${i}.oracle"
    fi
  done

  ## liblinear learning
  echo "liblinear learning ${root}weights.$iter" >&2
  qsubwrapper learn -l ${root}learn.$iter.log $cicada/cicada_learn_kbest \
                        --kbest  $tstset \
                        --oracle $orcset \
                        --output ${root}weights.$iter \
                        \
                        --learn-linear \
                        --solver 1 \
                        --C $C \
                        \
                        --debug=2 || exit 1


done 
