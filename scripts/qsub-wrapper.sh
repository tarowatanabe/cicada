#!/bin/sh
#
#  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#

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
hosts=""
hosts_file=""

### qsubs
mem_single=1gb
mem_mpi=8gb
queue=ltg

