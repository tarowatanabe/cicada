#!/bin/sh

me=`basename $0`

root=""
cicada=/data/lttools/decoder/cicada/bin
openmpi=""

tstset=""
refset=""

## # of processes, # of cores
np=4
nc=2
hosts=""
hosts_file=""

### decoding parameter
parse_cky=no

### MERT learning
iteration=10
weights_init=""
lower=""
upper=""

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
  Generation options
  --root                    root directory
  -c, --cicada              cicada directory (required)
  -m, --mpi                 MPI directory
  --host, --hosts           MPI hosts
  --hostfile, --host-file   MPI host file
  -q, --queue               PBS queue
  -n, --num                 # of processes to run

  Decoding options
  -p, --parse               perform parsing-composition with initial weights
  
  Training options
  -i, --iteration           MERT iterations
  -w, --weights             initial weights
  -l, --lower               lower-bound for features
  -u, --uppper              upper-bound for features
  -t, --test, --tstset      tuning data (required)
  -r, --reference, --refset reference translations (required)

  -h, --help                help message
"

while test $# -gt 0 ; do
  case $1 in
  --cicada | -c )
    test $# = 1 && eval "$exit_missing_arg"
    cicada=$2
    shift; shift ;;
  --mpi | -m )
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
  ## decoding
  --parse | -p )
    parse_cky=yes
    shift ;;

  ## training
  --iteration | -i )
    test $# = 1 && eval "$exit_missing_arg"
    iteration=$2
    shift; shift ;;
  --weights | -w )
    test $# = 1 && eval "$exit_missing_arg"
    weights_init=$2
    shift; shift ;;
  --lower | -l )
    test $# = 1 && eval "$exit_missing_arg"
    lower=$2
    shift; shift ;;
  --upper | -u )
    test $# = 1 && eval "$exit_missing_arg"
    lower=$2
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
  if test -e ${weights_init}; then
    weights="weights=${weights_init}"  
  fi

  if test -e "${root}weights.$iter_prev"; then
    weights="weights=${root}weights.$iter_prev"
  fi

  # generate translation and generate hypergraph (kbest=0)
  # instead of compose-cky, we use parse-cky with ${weights_init} learned by learn.sh!

  compose="compose-cky"
  if test "$parse_cky" = "yes"; then
    compose="parse-cky:size=1000,weights=${weights_init}"
  fi
 
  qsubwrapper decode -l ${root}decode.$iter.log ${openmpi}mpirun $mpinp $cicada/cicada_mpi \
	--input $tstset \
	--config ${root}cicada.config \
	\
	--operation $compose \
	--operation push-bos-eos \
	--operation apply:prune=true,size=200,${weights} \
	--operation prune:density=100,${weights} \
	--operation remove-bos-eos \
	--operation output:kbest=0,directory=${root}forest-$iter,${weights} \
	\
	--debug

  # compute bleu from hypergraph with current weights. this will also server as to how to grab 1best... see kbest=1 !
  qsubwrapper onebest -l ${root}1best.$iter.log ${openmpi}mpirun $mpinp $cicada/cicada_mpi \
	--input ${root}forest-$iter \
	--input-forest --input-directory \
	--operation output:kbest=1,${weights},file=${root}1best-$iter \
	--debug

  qsubsingle eval $cicada/cicada_eval --refset $refset --tstset ${root}1best-$iter --output ${root}eval-$iter.1best --scorer bleu:order=4

  ### forests upto now...
  tsts=""
  for ((i=1;i<=$iter;++i)); do
    if test -e ${root}forest-$i; then
      tsts="$tsts ${root}forest-$i"
    fi
  done

  ### previous weights...
  weights="${weights_init}"
  for ((i=1;i<$iter;++i)); do
    if test -e ${root}weights.$i; then
      weights="$weights ${root}weights.$i"
    fi
  done

  if test "$weights" != ""; then
    weights=" --feature-weights $weights"
  fi

  lower_bound=""
  if test "$lower" != ""; then
    lower_bound=" --bound-lower $lower"
  fi
  upper_bound=""
  if test "$upper" != ""; then
    upper_bound=" --bound-upper $upper"
  fi

  qsubwrapper mert -l ${root}mert.$iter.log ${openmpi}mpirun $mpinp $cicada/cicada_mert_mpi \
			--refset $refset \
			--tstset $tsts \
			--output ${root}weights.$iter \
			\
			$weights \
			--value-lower -5 \
			--value-upper  5 \
                        $lower_bound \
                        $upper_bound \
			\
			--normalize-l1 \
			--initial-average \
			\
			--debug=2

done 
