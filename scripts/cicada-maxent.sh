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
compose="compose-cky"

### linear learning
C=1e-3
oracle_cube=400

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
  -n, --np                  # of processes to run
  --nc                      # of cores to run
  
  Decoding options
  -c, --config              Configuration file (required)
  --compose                 Composition algorithm (default: compose-cky)
  
  Training options
  -C, --C                   hyperparameter
  --oracle-cube             cube size for oracle computation

  -d, --dev, --devset              tuning data (required)
  -r, --reference, --refset, --ref reference translations (required)

  -h, --help                help message
"

while test $# -gt 0 ; do
  case $1 in
  --root )
    test $# = 1 && eval "$exit_missing_arg"
    root=$2
    shift; shift ;;
  --cicada | --cicada-dir )
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
  --np | -n )
    test $# = 1 && eval "$exit_missing_arg"
    np=$2
    shift; shift ;;
  --nc )
    test $# = 1 && eval "$exit_missing_arg"
    nc=$2
    shift; shift ;;

  ## training
  --C | -C )
    test $# = 1 && eval "$exit_missing_arg"
    C=$2
    shift; shift ;;
  --oracle-cube )
    test $# = 1 && eval "$exit_missing_arg"
    oracle_cube=$2
    shift; shift ;;

  --config | -c )
    test $# = 1 && eval "$exit_missing_arg"
    config=$2
    shift; shift ;;
  --compose )
    test $# = 1 && eval "$exit_missing_arg"
    compose=$2
    shift; shift ;;

### test set and reference set
  --dev | -d | --devset )
    test $# = 1 && eval "$exit_missing_arg"
    devset=$2
    shift; shift ;;
  --reference | -r | --refset | --ref )
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
if test ! -e $devset; then
  echo "specify development data"
  exit 1
fi
if test "$refset" = ""; then
  echo "specify reference data"
  exit 1
fi
if test ! -e $refset; then
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

## check cicada...
cicadas="cicada_filter_config cicada_filter_weights cicada cicada_mpi cicada_eval cicada_oracle cicada_oracle_mpi cicada_learn cicada_learn_mpi"

found=yes
for prog in $cicadas; do
  if test ! -e "$cicada/$prog"; then
    found=no
    break
  fi
done

if test "$found" = "no"; then
  for bin in progs bin; do
    found=yes
    for prog in $cicadas; do
      if test ! -e "$cicada/$bin/$prog"; then
        found=no
        break
      fi
    done
    if test "$found" = "yes"; then
      cicada=$cicada/$bin
      break
    fi
  done
  
  if test "$found" = "no"; then
    echo "no --cicada | --cicada-dir?"
    exit 1
  fi
fi

if test "$openmpi" != ""; then
  openmpi=`echo "${openmpi}/" | sed -e 's/\/\/$/\//'`
fi

if test "$root" != ""; then
 root=`echo "${root}/" | sed -e 's/\/\/$/\//'`

  if test ! -e $root; then
    mkdir -p $root
  fi
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

  stripped=`expr "$1" : '\(.*\)_mpi$'`
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


### setup config file
### we will simply remove operation field..
echo "generate config file ${root}cicada.config.maxent" >&2
qsubwrapper config ${cicada}/cicada_filter_config \
      --remove-operation \
      --remove-feature-function \
      --input $config \
      --output ${root}cicada.config.maxent || exit 1

### actual composition
echo "composition ${root}forest-maxent" >&2
qsubwrapper decode -l ${root}forest.maxent.log $cicada/cicada_mpi \
	--input $devset \
	--config ${root}cicada.config.maxent \
        --operation $compose \
        --operation output:directory=${root}forest-maxent \
	\
	--debug || exit 1
  
echo "oracle translations ${root}forest-maxent.oracle" >&2
qsubwrapper oracle -l ${root}oracle.maxent.log $cicada/cicada_oracle_mpi \
        --refset $refset \
        --tstset ${root}forest-maxent \
        --output ${root}forest-maxent.oracle \
        --directory \
	--forest \
        --scorer  bleu:order=4,exact=true \
	--cube-size $oracle_cube \
        \
        --debug || exit 1

echo "learning ${root}weights.maxent" >&2
qsubwrapper learn -l ${root}learn.maxent.log $cicada/cicada_learn_mpi \
         --forest      ${root}forest-maxent \
         --intersected ${root}forest-maxent.oracle \
         --output ${root}weights.maxent \
         \
         --learn-lbfgs \
         --C $C \
         \
         --debug=2 || exit 1
