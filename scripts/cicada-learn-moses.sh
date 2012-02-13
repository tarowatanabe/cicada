#!/bin/sh
#
#  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#

### cicada-learn for moses...
### we will work with k-bests

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
moses=""
moses_options=""
moses_thread=no
config=""
### linear learning
iteration=20
iteration_first=1
weights_init=""
C=1e-3
regularize_l1=no
regularize_l2=no
scorer="bleu:order=4,exact=true"
learn="lbfgs"
learn_options=""
zero_weights=no
history_weights=no
kbest=1000
bias_features=""
bias_weight=-1
merge="no"
interpolate=0.0
lower=""
upper=""

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
  --moses                   moses-cmd binary
  --mpi                     MPI directory
  --host, --hosts           MPI hosts
  --hostfile, --host-file   MPI host file
  -q, --queue               PBS queue                (default: $queue)
  -n, --np                  # of processes to run    (default: $np)
  --nc                      # of cores to run        (default: $nc)
  --mem                     memory used by each node (default: $mem)

  Decoding options
  -c, --config              Configuration file (required)
  -o, --options             Moses options
  --thread                  run moses with multiple threads
  
  Training options
  -i, --iteration           PRO iterations    (default: $iteration)
  --iteration-first         The first iteration (default: $iteration_first)
  -w, --weights             initial weights
  -C, --C                   hyperparameter    (default: $C)
  --regularize-l1           L1 regularization
  --regularize-l2           L2 regularization                (default)
  --scorer                  scorer            (default: $scorer)
  --learn                   learner (lbfgs, svm, linear, sgd, pegasos, mira, cw, arow, nherd, cp, mcp, xbleu)
                            (WARNING: --learn-liner or --liblinear option is deprecated. use --learn linear)
  --learn-options           other learning options
  --zero-weights            learning from zero weights in each iteration
  --history-weights         learning from multiple hisories (only for mcp learning)
  --kbest                   kbest size             (default: $kbest)
  --bias-features           bias features
  --bias-weight             bias weight
  --merge                   perform kbest merging
  --interpolate             weights interpolation

  --lower                   lower-bound for features
  --uppper                  upper-bound for features

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
  --moses )
    test $# = 1 && eval "$exit_missing_arg"
    moses=$2
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
  --iteration-first )
    test $# = 1 && eval "$exit_missing_arg"
    iteration_first=$2
    shift; shift ;;
  --weights | -w )
    test $# = 1 && eval "$exit_missing_arg"
    weights_init=$2
    shift; shift ;;
  --C | -C )
    test $# = 1 && eval "$exit_missing_arg"
    C=$2
    shift; shift ;;
  --regularize-l1 )
    regularize_l1=yes
    shift ;;
  --regularize-l2 )
    regularize_l2=yes
    shift ;;
  --scorer )
    test $# = 1 && eval "$exit_missing_arg"
    scorer=$2
    shift; shift ;;
  --learn )
    test $# = 1 && eval "$exit_missing_arg"
    learn=$2
    shift; shift ;;
  --learn-options )
    test $# = 1 && eval "$exit_missing_arg"
    learn_options=$2
    shift; shift ;;
  --zero-weights )
    zero_weights=yes
    shift ;;
  --history-weights )
    history_weights=yes
    shift ;;

  --kbest )
    test $# = 1 && eval "$exit_missing_arg"
    kbest=$2
    shift; shift ;;
  --bias-features )
    test $# = 1 && eval "$exit_missing_arg"
    bias_features=$2
    shift; shift ;;
  --bias-weight )
    test $# = 1 && eval "$exit_missing_arg"
    bias_weight=$2
    shift; shift ;;

  --merge )
    merge=yes
    shift ;;
  --interpolate )
    test $# = 1 && eval "$exit_missing_arg"
    interpolate=$2
    shift; shift ;;

  --lower )
    test $# = 1 && eval "$exit_missing_arg"
    lower=$2
    shift; shift ;;
  --upper )
    test $# = 1 && eval "$exit_missing_arg"
    lower=$2
    shift; shift ;;

  --config | -c )
    test $# = 1 && eval "$exit_missing_arg"
    config=$2
    shift; shift ;;
  --options | -o )
    test $# = 1 && eval "$exit_missing_arg"
    moses_options=$2
    shift; shift ;;
  --thread )
    moses_thread=yes
    shift ;;

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

if test "$devset" = "" -o ! -e "$devset"; then
  echo "specify development data" >&2
  exit 1
fi
if test "$refset" = "" -o ! -e "$refset"; then
  echo "specify reference data" >&2
  exit 1
fi
if test "$config" = "" -o ! -e "$config"; then
  echo "specify config file" >&2
  exit 1
fi
if test "$cicada" = ""; then
  echo "no cicada dir?" >&2
  exit 1
fi
if test "$moses" = "" -o ! -x "$moses"; then
  echo "no moses" >&2
  exit 1
fi

if test "$regularize_l1" = no -a "$regularize_l2" = no; then
  regularize_l2=yes
fi

if test "$regularize_l1" = yes -a "$regularize_l2" = yes; then
  echo "both L1 and L2?" >&2
  exit 1  
fi

learner="cicada_learn_kbest_mpi"
learn_option=""
case $learn in
  lbfgs | sgd | pegasos | mira | cw | arow | nherd | cp | mcp | xbleu )
    learner="cicada_learn_kbest_mpi"
    learn_option=" --learn-$learn"
    break ;;
  linear )
    learner="cicada_learn_kbest"
    learn_option=" --learn-linear"
    break ;;
  svm )
    learner="cicada_learn_kbest"
    learn_option=" --learn-svm"
    break ;;
  * )
    echo "learning algorithm can be either lbfgs, linear, svm, sgd, pegasos, mira, cw, arow, cp, mcp"
    exit 1 ;;
esac

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

## check cicada...
cicadas="cicada_filter_config_moses cicada_filter_kbest_moses cicada_filter_weights cicada_eval cicada_oracle_kbest cicada_learn_kbest cicada-moses.sh"

for prog in $cicadas; do
  tmptmp=`cicadapath $prog`
  
  if test ! -e $tmptmp; then
    echo "no $prog at $tmptmp... no --cicada | --cicada-dir?" >&2
    exit 1
  fi
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
  args__=""
  for arg in "$@"; do 
    args__="$args__ `argument $arg`"
  done
  echo $args__
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
      
      if test "$mpimode" = "yes"; then
	echo "${openmpi}mpirun $mpinp `arguments "$@"` $out_option $log_option"
      else
	## shift here!
	shift;
	echo "$stripped `arguments "$@"` $threads $out_option $log_option"
      fi
    ) |
    qsub -S /bin/sh || exit 1
  else
    if test "$mpimode" = "yes"; then
      eval "${openmpi}mpirun $mpinp `arguments "$@"` $out_option $log_option" || exit 1
    else
      shift
      eval "$stripped `arguments "$@"` $threads $out_option $log_option" || exit 1
    fi
  fi
}


for ((iter=$iteration_first;iter<=iteration; ++ iter)); do
  echo "iteration: $iter" >&2
  iter_prev=`expr $iter - 1`
  
  ## START of moses specific changes...
  
  ### generate moses.ini file, using $config and weights.$iter_prev
  ### HOW?
  
  moses_ini=${root}moses.ini.$iter
  
  weights_process=""
  if test "${weights_init}" != ""; then
    weights_process=$weights_init
  fi
  
  if test -e ${root}weights.$iter_prev; then
    weights_process=${root}weights.$iter_prev
  fi

  if test "$weights_process" = ""; then
    qsubwrapper weights `cicadapath cicada_filter_config_moses` \
	--input $config \
        --output weights.moses.init
    
    weights_process=weights.moses.init
  fi
    
  if test "$bias_features" != ""; then
    qsubwrapper config `cicadapath cicada_filter_config_moses` \
	--weights $weights_process \
	--bias-features $bias_features \
	--bias-weight   $bias_weight \
	--input $config \
	--output $moses_ini || exit 1
  else
    qsubwrapper config `cicadapath cicada_filter_config_moses` \
	--weights $weights_process \
	--input $config \
	--output $moses_ini || exit 1
  fi
  
  ### run cicada wrapped moses...
  output=${root}kbest-$iter

  ### create output directory
  if test -e $output; then
    rm -rf $output || exit 1
  fi
  if test ! -e $output; then
    mkdir -p $output || exit 1
  fi

  echo "decoding $output" >&2

  if test "$moses_thread" = yes; then
    
    kbest_file=${root}kbest.$kbest.$iter
    
    qsubwrapper moses -l ${root}decode.$iter.log -o ${root}decode.$iter.out \
	$moses \
	-input-file $devset \
	-config $moses_ini \
	$moses_options \
	-n-best-list $kbest_file $kbest distinct \
	-threads $nc || exit 1
    
    # kbest filtering..
    if test "$bias_features" != ""; then
      qsubwrapper kbest \
	`cicadapath cicada_filter_kbest_moses` \
	--input $kbest_file \
	--output $output \
        --erase-features $bias_features \
	--directory || exit 1
    else
      qsubwrapper kbest \
	`cicadapath cicada_filter_kbest_moses` \
	--input $kbest_file \
	--output $output \
	--directory || exit 1
    fi
  else
    mkdir -p $output/kbests || exit 1

    ### generate scripts for kbest generation
    kbest_generation=${output}/kbests/kbest-generation
    for ((i=0;i<$np;++i)); do
      kbest_file=${output}/kbests/kbest.$i

      filter=`cicadapath cicada_filter_kbest_moses`
  
      kbest_option="-n-best-list $kbest_file $kbest distinct"
      moses_cmd="$moses -config $moses_ini $moses_options $kbest_option"
      
      filter_cmd="$filter --input $kbest_file --output $output --directory --keep --offset $i --stride $np "
      if test "$bias_features" != ""; then
        filter_cmd="$filter_cmd --erase-features $bias_features"
      fi
    
      echo  "$moses_cmd && $filter_cmd" >> $kbest_generation
    done  
  
    ### actually run
    qsubwrapper kbest -m -l ${root}decode.$iter.log -o ${root}decode.$iter.out \
        `cicadapath mpimap` \
        --prog `cicadapath mpimap` \
        --even \
        --input $devset \
        $kbest_generation || exit 1
    
    ## remove unused directory...
    rm -rf ${output}/kbests || exit 1
  fi

  ### END of moses specific changes...

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

  ### previous weights...
  weights_last=${weights_init}
  weights_history=${weights_init}
  for ((i=1;i<$iter;++i)); do
    if test -e ${root}weights.$i; then
      weights_last=${root}weights.$i
      weights_history="${weights_history} ${weights_last}"
    fi
  done

  ### option for previous weights
  weights_option=""
  if test "$weights_last" != ""; then
    weights_option=" --weights $weights_last"
  fi
  if test "$zero_weights" = "yes"; then
    weights_option=""
  fi
  if test "$history_weights" = "yes"; then
    if test "$weights_history" != ""; then
      weights_option="${weights_option} --weights-history ${weights_history}"
    fi
  fi

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

  regularize=" --regularize-l2"
  if test "$regularize_l1" = yes; then
    regularize=" --regularize-l1"
  fi

  lower_bound=""
  if test "$lower" != ""; then
    lower_bound=" --bound-lower $lower"
  fi
  upper_bound=""
  if test "$upper" != ""; then
    upper_bound=" --bound-upper $upper"
  fi

  echo "learning ${root}weights.$iter" >&2
  qsubwrapper learn -t -l ${root}learn.$iter.log `cicadapath $learner` \
                        --kbest  $tstset \
                        --oracle $learn_oracle \
                        --refset $refset \
                        --scorer $scorer \
	                $unite \
                        --output $weights_learn \
                        \
                        $weights_option \
                        $learn_option \
                        $learn_options \
                        $regularize \
                        --C $C \
                        $lower_bound \
                        $upper_bound \
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
