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

### linear learning
iteration=20
weights_init=""
C=1e-3
scorer="bleu:order=4,exact=true"
solver=1
kbest=1000
merge="no"
interpolate=0.0

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
  
  Training options
  -i, --iteration           PRO iterations    (default: $iteration)
  -w, --weights             initial weights
  -C, --C                   hyperparameter    (default: $C)
  --scorer                  scorer            (default: $scorer)
  --solver                  liblinear solver type. See liblinear FAQ,
                            or run cicada_learn_kbest --help
                            (Default: 1, L2-reg, L2-loss SVM)
  --kbest                   kbest size             (default: $kbest)
  --merte                   perform kbest merging
  --interpolate             weights interpolation

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
  --scorer )
    test $# = 1 && eval "$exit_missing_arg"
    scorer=$2
    shift; shift ;;
  --solver )
    test $# = 1 && eval "$exit_missing_arg"
    solver=$2
    shift; shift ;;
  --kbest )
    test $# = 1 && eval "$exit_missing_arg"
    kbest=$2
    shift; shift ;;
  --merge )
    merge=yes
    shift ;;
  --interpolate )
    test $# = 1 && eval "$exit_missing_arg"
    interpolate=$2
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

cicadas="cicada_filter_config cicada_filter_weights cicada cicada_mpi cicada_eval cicada_oracle_kbest cicada_oracle_kbest_mpi cicada_learn_kbest"

for prog in $cicadas; do
  tmptmp=`cicadapath $prog` || (echo "no $prog... no --cicada | --cicada-dir?" >&2; exit 1)
done

if test "$weights_init" != ""; then
  if test ! -e $weights_init; then
    echo "no initial weights: $weights_init ?" >&2
    exit 1
  fi
fi

if test "$openmpi" != ""; then
  openmpi=`echo "${openmpi}/" | sed -e 's/\/\/$/\//'`

  if test ! -e ${openmpi}mpirun; then
    openmpi=${openmpi}bin/
    if test ! -e ${openmpi}mpirun; then
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

do_interpolate=`echo "($interpolate > 0.0) && ($interpolate < 1.0)" | bc`

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

argument() {
  if test $# -gt 1; then
    echo "\"$@\""
  else
    echo "$@"
  fi
}

arguments() {
  args=""
  for arg in "$@"; do 
    args="$args `argument $arg`"
  done
  echo $args
}

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

  out_option=""
  if test "$outfile" != ""; then
    out_option="> $outfile"
  fi
  log_option=""
  if test "$logfile" != ""; then
    log_option="2> $logfile"
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
      
      ### we need to handle argument spiltting...
      if test "$mpimode" = "yes"; then
        parameters=`arguments "$@"`
	echo "${openmpi}mpirun $mpinp $parameters $out_option $log_option"
      else
	## shift here!
	shift;
	parameters=`arguments "$@"`
	echo "$stripped $parameters $threads $out_option $log_option"
      fi
    ) |
    qsub -S /bin/sh || exit 1
  else
    if test "$mpimode" = "yes"; then
      parameters=`arguments "$@"`
      eval "${openmpi}mpirun $mpinp $parameters $out_option $log_option" || exit 1
    else
      shift
      parameters=`arguments "$@"`
      eval "$stripped $parameters $threads $out_option $log_option" || exit 1
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
  qsubwrapper config `cicadapath cicada_filter_config` \
      --weights $weights \
      --kbest   $kbest \
      --file directory=${root}kbest-$iter \
      --input $config \
      --output ${root}cicada.config.$iter || exit 1
  
  ### actual decoding
  echo "decoding ${root}kbest-$iter" >&2
  qsubwrapper decode -l ${root}decode.$iter.log `cicadapath cicada_mpi` \
	--input $devset \
	--config ${root}cicada.config.$iter \
	\
	--debug || exit 1

  ### BLEU
  echo "BLEU ${root}eval-$iter.1best" >&2
  qsubwrapper eval `cicadapath cicada_eval` \
      --refset $refset \
      --tstset ${root}kbest-$iter \
      --output ${root}eval-$iter.1best \
      --scorer $scorer || exit 1

  ### kbests upto now...
  tstset=""
  orcset=""
  for ((i=1;i<=$iter;++i)); do
    if test -e ${root}kbest-$i; then
      tstset="$tstset ${root}kbest-$i"
      orcset="$orcset ${root}kbest-${i}.oracle"
    fi
  done

  ### last weights
  weights_last=${weights_init}
  for ((i=1;i<$iter;++i)); do
    if test -e ${root}weights.$i; then
      weights_last=${root}weights.$i
    fi
  done

  tstset_oracle=${root}kbest-$iter
  if test "$merge" = "yes"; then
    tstset_oracle=$tstset
  fi

  ### compute oracles
  echo "oracle translations ${root}kbest-${iter}.oracle" >&2
  qsubwrapper oracle -t -l ${root}oracle.$iter.log `cicadapath cicada_oracle_kbest_mpi` \
        --refset $refset \
        --tstset $tstset_oracle \
        --output ${root}kbest-${iter}.oracle \
        --directory \
        --scorer  $scorer \
        \
        --debug || exit 1

  learn_oracle=$orcset
  unite=""
  if test "$merge" = "yes"; then
    learn_oracle=${root}kbest-${iter}.oracle
    unite=" --unite"
  fi

  weights_learn=${root}weights.$iter
  if test $do_interpolate -eq 1; then
    weights_learn=${root}weights.${iter}.learn
  fi
    
  ## liblinear learning
  echo "liblinear learning ${root}weights.$iter" >&2
  qsubwrapper learn -t -l ${root}learn.$iter.log `cicadapath cicada_learn_kbest` \
                        --kbest  $tstset \
                        --oracle $learn_oracle \
	                $unite \
                        --output $weights_learn \
                        \
                        --learn-linear \
                        --solver $solver \
                        --C $C \
                        \
                        --debug=2 || exit 1

  ### interpolate from the previous iterations, if interpolate_ratio is >0 and <1
  if test $do_interpolate -eq 1; then
    if test "$weights_last" = ""; then
      cp $weights_learn ${root}weights.$iter
    else
      echo "merging weights" >&2
      qsubwrapper interpolate `cicadapath cicada_filter_weights` \
	  --output ${root}weights.$iter \
	  ${weights_last}:scale=`echo "1.0 - $interpolate_ratio" | bc`  \
          ${weights_learn}:scale=$interpolate_ratio || exit 1
    fi
  fi

done 
