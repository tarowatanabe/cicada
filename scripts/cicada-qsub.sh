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

## # of processes, # of cores
np=1
nc=2
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
  -q, --queue               PBS queue
  -n, --np                  # of processes to run
  --nc                      # of cores to run
  --mem                     memory used by each node (default: $mem)
  --log, -l                 logfile
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
  --log | -l )
    test $# = 1 && eval "$exit_missing_arg"
    logfile=$2
    shift; shift ;;
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

if test "$openmpi" != ""; then
  openmpi=`echo "${openmpi}/" | sed -e 's/\/\/$/\//'`
fi


if test "$qsub" != ""; then
    
  if test "$name" = ""; then
    echo "no process name under qsub (specify by --name)"
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
      echo "#PBS -l select=$np:ncpus=$nc:mpiprocs=$nc:mem=${mem}"
      echo "#PBS -l place=scatter"
    else
      echo "#PBS -l select=1:ncpus=1:mem=${mem}"
    fi
      
    if test "$TMPDIR_SPEC" != ""; then
      echo "export TMPDIR_SPEC=$TMPDIR_SPEC"
    fi
    if test "$TMPDIR" != ""; then
      echo "export TMPDIR=$TMPDIR"
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
  exec qsub -S /bin/sh
else
  if test "$stripped" != "$1" -a $np -gt 1; then
    if test "$logfile" != ""; then
      exec ${openmpi}mpirun $mpinp $@ >& $logfile
    else
      exec ${openmpi}mpirun $mpinp $@
    fi
  else
    shift
    if test "$logfile" != ""; then
      exec $stripped $@ >& $logfile
    else
      exec $stripped $@
    fi
  fi
fi
