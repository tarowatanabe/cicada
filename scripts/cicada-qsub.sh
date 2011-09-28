#!/bin/sh
#
#  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#

### we assume PBSPro. If you want to apply this to other environgmnet, adjust 
### #PBS stuff and qsub related commands

me=`basename $0`

exit_missing_arg="\
echo \"$me: option \\\`\$1' requires an argument\" >&2
echo \"\$help\" >&2
exit 1"

### working dir..
workingdir=`pwd`

### this is a test, whether we can run under cluster or not...
qsub=`which qsub 2> /dev/null`

name=""
logfile=""
outfile=""
threads=""

## # of processes, # of cores
np=1
nc=1
openmpi=""
hosts=""
hosts_file=""

### qsubs
mem=1gb
queue=ltg

usage="\
$me name [options] cicada-program args
  --mpi                     MPI implementation
  --host, --hosts           MPI hosts
  --hostfile, --host-file   MPI host file
  -q, --queue               PBS queue                (default: $queue)
  -n, --np                  # of processes to run    (default: $np)
  --nc                      # of cores to run        (default: $nc)
  --mem                     memory used by each node (default: $mem)
  --log, --logfile, -l      logfile
  --out, --outfile, -o      stdout file
  --threads,-t              threading option
  --name                    process name

  -h, --help                help message
"

while test $# -gt 0; do
  case $1 in
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
  --log | --logfile | -l )
    test $# = 1 && eval "$exit_missing_arg"
    logfile=$2
    shift; shift ;;
  --out | --outfile | -o )
    test $# = 1 && eval "$exit_missing_arg"
    outfile=$2
    shift; shift ;;
  --threads | -t )
    threads="yes"
    shift ;;
  --name )
    test $# = 1 && eval "$exit_missing_arg"
    name=$2
    shift; shift ;;

  --help | -h )
    echo "$usage"
    exit ;;
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

if test "$threads" = "yes"; then
  threads=" --threads ${nc}"
fi

if test "$openmpi" != ""; then
  openmpi=`echo "${openmpi}/" | sed -e 's/\/\/$/\//'`
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


if test "$qsub" != ""; then
    
  if test "$name" = ""; then
    echo "no process name under qsub (specify by --name)" >&2
    exit 1
  fi

  (
    echo "#!/bin/sh"
    echo "#PBS -N $name"
    echo "#PBS -W block=true"
    echo "#PBS -e /dev/null"
    echo "#PBS -o /dev/null"
    echo "#PBS -q $queue"
    if test "$stripped" != "$1" -a $np -gt 1; then
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

    if test "$stripped" != "$1" -a $np -gt 1; then
      if test "$logfile" != ""; then
	if test "$outfile" != ""; then
          echo "${openmpi}mpirun $mpinp $@ > $outfile 2> $logfile"
	else
	  echo "${openmpi}mpirun $mpinp $@ 2> $logfile"
	fi
      else
	if test "$outfile" != ""; then
          echo "${openmpi}mpirun $mpinp $@ > $outfile"
	else
	  echo "${openmpi}mpirun $mpinp $@"
	fi
      fi
    else
      ## shift here!
      shift;
      if test "$logfile" != ""; then
	if test "$outfile" != ""; then
          echo "$stripped $@ $threads > $outfile 2> $logfile"
	else
	  echo "$stripped $@ $threads 2> $logfile"
	fi
      else
	 if test "$outfile" != ""; then
          echo "$stripped $@ $threads > $outfile"
	else
	  echo "$stripped $@ $threads"
	fi
      fi
    fi
  ) |
  exec qsub -S /bin/sh
else
  if test "$stripped" != "$1" -a $np -gt 1; then
    if test "$logfile" != ""; then
      if test "$outfile" != ""; then
        exec ${openmpi}mpirun $mpinp "$@" > $outfile 2> $logfile
      else
        exec ${openmpi}mpirun $mpinp "$@" 2> $logfile
      fi
    else
      if test "$outfile" != ""; then
        exec ${openmpi}mpirun $mpinp "$@" > $outfile
      else
	exec ${openmpi}mpirun $mpinp "$@"
      fi
    fi
  else
    shift
    if test "$logfile" != ""; then
      if test "$outfile" != ""; then
	exec $stripped "$@" $threads > $outfile 2> $logfile
      else
	exec $stripped "$@" $threads 2> $logfile
      fi
    else
      if test "$outfile" != ""; then
        exec $stripped "$@" $threads > $outfile
      else
	exec $stripped "$@" $threads
      fi
    fi
  fi
fi
