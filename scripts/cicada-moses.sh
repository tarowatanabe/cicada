#!/bin/sh
#
#  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#

#
# cicada wrapper for moses :-)
#
# we support multiple-node execution
#

### single-best mode
### use mpipe for grabbing stdout from moses
### 

### 
### kbest mode (we need additional parameters, # of procs...)
###
### we use mpimap for mapping input evenly (map), and reduce results by filter_kbest_moses (and delete moses-output)
###
### first, directory is created
### second, generate moses-shell-command into the direcotry,
### third, run mpimap with moses-shell-commands
### forth, collect moses outputs
### fifth, erase unused files!
###

me=`basename $0`

### working dir..
workingdir=`pwd`

### this is a test, whether we can run under cluster or not...
qsub=`which qsub 2> /dev/null`

### np, nc config
cicada=""
openmpi=""

## # of processes, # of cores
np=1
nc=1
hosts=""
hosts_file=""
mem=8gb
queue=ltg

### decoding config
input="-"
output="-"
moses=""
moses_config=""
moses_options=""

### kbest options
kbest=1
kbest_unique="no"


exit_missing_arg="\
echo \"$me: option \\\`\$1' requires an argument\" >&2
echo \"\$help\" >&2
exit 1"

usage="\
$me [options]
  General options
  --cicada                  cicada directory (required)
  --moses                   moses-cmd binary
  --mpi                     MPI directory
  --host, --hosts           MPI hosts
  --hostfile, --host-file   MPI host file
  -q, --queue               PBS queue                (default: $queue)
  -n, --np                  # of processes to run    (default: $np)
  --nc                      # of cores to run        (default: $nc)
  --mem                     memory used by each node (default: $mem)

  Input/Output
  --input                   input file  (default: $input)
  --output                  output file (default: $output)

  Decoding options
  -c, --config              Configuration file (required)
  -o, --options             Moses options
  --kbest                   kbest size             (default: $kbest, meaning single best)
  --unique                  unique kbest

  -h, --help                help message
"

while test $# -gt 0 ; do
  case $1 in
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

  --input )
    test $# = 1 && eval "$exit_missing_arg"
    input=$2
    shift; shift ;;
  --output )
    test $# = 1 && eval "$exit_missing_arg"
    output=$2
    shift; shift ;;

  --config | -c )
    test $# = 1 && eval "$exit_missing_arg"
    moses_config=$2
    shift; shift ;;
  --options | -o )
    test $# = 1 && eval "$exit_missing_arg"
    moses_options=$2
    shift; shift ;;

  --kbest )
    test $# = 1 && eval "$exit_missing_arg"
    kbest=$2
    shift; shift ;;
  --unique )
    kbest_unique=yes
    shift ;;

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

if test "$moses_config" = ""; then
  echo "specify config file" >&2
  exit 1
fi
if test "$cicada" = ""; then
  echo "no cicada dir?" >&2
  exit 1
fi

if test ! -e $moses_config; then
  echo "no config file?" >& 2
  exit 1
fi
if test ! -x $moses; then
  echo "no moses" >&2
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

cicadas="cicada_filter_kbest_moses mpimap mpipe mpish"

for prog in $cicadas; do
  tmptmp=`cicadapath $prog` || (echo "no $prog... no --cicada | --cicada-dir?" >&2; exit 1)
done

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

if test $kbest -le 1; then
  ### 1best mode
  ### we will use mpipe
  ### TODO: this will not work for PBS, since we need to evaluate this again...
    
  command="$moses -config $moses_config $moses_options"
  
  qsubwrapper 1best -m `cicadapath mpipe` --command "$command" --input $input --output $output --debug || exit 1
  
else
  ### kbest mode

  if test "$output" = "-" -o "$output" = ""; then
    echo "no output file" >&2
    exit 1
  fi
    
  ### create output directory
  if test -e $output; then
    rm -rf $output || exit 1
  fi
  if test ! -e $output; then
    mkdir -p $output || exit 1
  fi
  
  mkdir -p $output/kbests || exit 1
  
  ### generate scripts for kbest generation
  kbest_generation=${output}/kbests/kbest-generation
  for ((i=0;i<$np;++i)); do
    kbest_file=${output}/kbests/kbest.$i
  
    kbest_option="-n-best-list $kbest_file $kbest"
    if test "$kbest_unique" = "yes"; then
      kbest_option="$kbest_option distinct"
    fi

    filter=`cicadapath cicada_filter_kbest_moses`
 
    moses_cmd="$moses -config $moses_config $moses_options $kbest_option"
    filter_cmd="$filter --input $kbest_file --output $output --directory --keep --offset $i --stride $np"
    
    echo "$moses_cmd && $filter_cmd" >> $kbest_generation
  done
  
  ### actually run
  qsubwrapper kbest -m `cicadapath mpimap` \
      --prog `cicadapath mpimap` \
      --even \
      --input $input \
      $kbest_generation || exit 1
  
  rm -rf ${output}/kbests || exit 1
fi
