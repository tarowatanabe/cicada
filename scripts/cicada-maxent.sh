#!/bin/sh
#
#  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#

### we assume PBSPro. If you want to apply this to other environgmnet, adjust 
### #PBS stuff and qsub related commands

me=`basename $0`

### working dir..
workingdir=`pwd`

### this is a test, whether we can run under cluster or not...
qsub=`which qsub 2> /dev/null`

root=""
cicada=""
openmpi=""

devset=""
refset=""

## # of processes, # of cores
np=1
nc=1
hosts=""
hosts_file=""

### decoding config
config=""
compose="compose-cky"

### linear learning
C=1e-3
oracle_cube=400
scorer="bleu:order=4,exact=true"

### qsubs
mem=8gb
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
  -q, --queue               PBS queue                (default: $queue)
  -n, --np                  # of processes to run    (default: $np)
  --nc                      # of cores to run        (default: $nc)
  --mem                     memory used by each node (default: $mem)
  
  Decoding options
  -c, --config              Configuration file (required)
  --compose                 Composition algorithm (default: $compose)
  
  Training options
  -C, --C                   hyperparameter                   (default: $C)
  --oracle-cube             cube size for oracle computation (default: $oracle_cube)
  --scorer                  scorer                           (default: $scorer)

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
  --mpi | --mpi-dir )
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
  --mem )
    test $# = 1 && eval "$exit_missing_arg"
    mem=$2
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
  --scorer )
    test $# = 1 && eval "$exit_missing_arg"
    scorer=$2
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
  echo "specify development data" >&2
  exit 1
fi
if test ! -e $devset; then
  echo "specify development data" >&2
  exit 1
fi
if test "$refset" = ""; then
  echo "specify reference data" >&2
  exit 1
fi
if test ! -e $refset; then
  echo "specify reference data" >&2
  exit 1
fi

if test "$config" = ""; then
  echo "specify config file" >&2
  exit 1
fi
if test "$cicada" = ""; then
  echo "no cicada dir?" >&2
  exit 1
fi

## check cicada...
cicadapath() {
  file=$1
  shift
  
  path=$cicada/$file
  if test ! -e $path; then
    path=$cicada/bin/$file
    if test ! -e $path; then
      path=$cicada/progs/$file
      if test ! -e $path; then
        path=$cicada/scripts/$file
	if test ! -e $path; then
	  echo $file
	  return 1
	fi
      fi
    fi
  fi
  echo $path
  return 0
}

cicadas="cicada_filter_config cicada_filter_weights cicada cicada_mpi cicada_eval cicada_oracle cicada_oracle_mpi cicada_learn cicada_learn_mpi"

for prog in $cicadas; do
  tmptmp=`cicadapath $prog` || (echo "no $prog... no --cicada | --cicada-dir?" >&2; exit 1)
done

if test "$openmpi" != ""; then
  openmpi=`echo "${openmpi}/" | sed -e 's/\/\/$/\//'`
  
  if test ! -e $openmpi/mpirun; then
    openmpi=$openmpi/bin
    if test ! -e $openmpi/mpirun; then
      echo "no mpirun?" >&2
    exit 1
    fi
  fi
fi

if test "$root" != ""; then
 root=`echo "${root}/" | sed -e 's/\/\/$/\//'`

  if test ! -e $root; then
    mkdir -p $root
  fi
fi

### check np and nc
if test $np -le 1; then
  np=1
fi
if test $nc -le 1; then
  nc=1
fi

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
  outfile=""
  threads=""
  mpimode=no
  while test $# -gt 0 ; do
  case $1 in
  -t )
    threads=" --threads ${nc}"
    shift ;;
  -m )
    mpimode=yes
    shift ;;
  -l )
    test $# = 1 && eval "$exit_missing_arg"
    logfile=$2
    shift; shift ;;
  -o )
    test $# = 1 && eval "$exit_missing_arg"
    outfile=$2
    shift; shift ;;
  -* )
    exec >&2
    echo "$me: invalid option $1" >&2
    exit 1 ;;
  * )
    break ;;
  esac
  done

  stripped=`expr "$1" : '\(.*\)_mpi$'`
  if test "$stripped" = ""; then
    stripped=$1
  fi

  if test "$mpimode" = "no"; then
    if test "$stripped" != "$1" -a $np -gt 1; then
      mpimode=yes
    fi
  fi

  if test "$qsub" != ""; then
    (
      echo "#!/bin/sh"
      echo "#PBS -N $name"
      echo "#PBS -W block=true"
      echo "#PBS -e /dev/null"
      echo "#PBS -o /dev/null"
      echo "#PBS -q $queue"
      if test "$mpimode" = "yes"; then
        echo "#PBS -l select=${np}:ncpus=${nc}:mpiprocs=${nc}:mem=${mem}"
        echo "#PBS -l place=scatter"
      else
        echo "#PBS -l select=1:ncpus=${nc}:mem=${mem}"
      fi

      if test "$TMPDIR_SPEC" != ""; then
        echo "export TMPDIR_SPEC=$TMPDIR_SPEC"
      fi
      if test "$TMPDIR" != ""; then
        echo "export TMPDIR=$TMPDIR"
      fi
      if test "$LD_LIBRARY_PATH" != ""; then
        echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
      fi

      echo "cd $workingdir"

      out_option=""
      if test "$outfile" != ""; then
	out_option="> $outfile"
      fi
      log_option=""
      if test "$logfile" != ""; then
	log_option="2> $logfile"
      fi
      
      ### we need to handle argument spiltting...
      if test "$mpimode" = "yes"; then
	echo "${openmpi}mpirun $mpinp $@ $out_option $log_option"
      else
	## shift here!
	shift;
	echo "$stripped $@ $threads $out_option $log_option"
      fi
    ) |
    qsub -S /bin/sh || exit 1
  else
    
    if test "$logfile" = ""; then
      logfile=/dev/stderr
    fi
    if test "$outfile" = ""; then
      outfile=/dev/stdout
    fi

    if test "$mpimode" = "yes"; then
      ${openmpi}mpirun $mpinp "$@" > $outfile 2> $logfile || exit 1
    else
      shift
      $stripped "$@" $threads > $outfile 2> $logfile || exit 1
    fi
  fi
}

### setup config file
### we will simply remove operation field..
echo "generate config file ${root}cicada.config.maxent" >&2
qsubwrapper config `cicadapath cicada_filter_config` \
      --remove-operation \
      --remove-feature-function \
      --input $config \
      --output ${root}cicada.config.maxent || exit 1

### actual composition
echo "composition ${root}forest-maxent" >&2
qsubwrapper decode -l ${root}forest.maxent.log `cicadapath cicada_mpi` \
	--input $devset \
	--config ${root}cicada.config.maxent \
        --operation $compose \
        --operation output:directory=${root}forest-maxent \
 	\
	--debug || exit 1
  
echo "oracle translations ${root}forest-maxent.oracle" >&2
qsubwrapper oracle -t -l ${root}oracle.maxent.log `cicadapath cicada_oracle_mpi` \
        --refset $refset \
        --tstset ${root}forest-maxent \
        --output ${root}forest-maxent.oracle \
        --directory \
	--forest \
        --scorer  $scorer \
	--cube-size $oracle_cube \
        \
        --debug || exit 1

echo "learning ${root}weights.maxent" >&2
qsubwrapper learn -t -l ${root}learn.maxent.log `cicadapath cicada_learn_mpi` \
         --forest ${root}forest-maxent \
         --oracle ${root}forest-maxent.oracle \
         --output ${root}weights.maxent \
         \
         --learn-lbfgs \
         --C $C \
         \
         --debug=2 || exit 1
