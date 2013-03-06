#!/usr/bin/env python
#
#  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
#
###
### a wrapper script for running multiple commands
### inspired by mpipe and thrpe, we support pbs 
### Actually, we will run by mpipe and thrpe

import threading
import multiprocessing

import time
import sys
import os, os.path
import string
import re
import subprocess
import shutil

import UserList
import UserString
import cStringIO

from optparse import OptionParser, make_option

opt_parser = OptionParser(
    option_list=[
    make_option("--input", default="-", action="store", type="string",
                metavar="FILE", help="input file"),
    make_option("--output", default="-", action="store", type="string",
                metavar="FILE", help="output file"),
    make_option("--logfile", default="", action="store", type="string",
                metavar="FILE", help="log file for stderr"),
    
    ### actual command to run!
    make_option("--command", default="", action="store", type="string",
                metavar="COMMAND", help="command to run"),
    
    ## max-malloc
    make_option("--max-malloc", default=4, action="store", type="float",
                metavar="MALLOC", help="maximum memory in GB (default: 4)"),

    # CICADA Toolkit directory
    make_option("--cicada-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="cicada directory"),
    # MPI Implementation.. if different from standard location...
    make_option("--mpi-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="MPI directory"),

    make_option("--threads", default=1, action="store", type="int",
                help="# of thrads for thread-based parallel processing"),
    # perform threading or MPI training    
    make_option("--mpi", default=0, action="store", type="int",
                help="# of processes for MPI-based parallel processing. Identical to --np for mpirun"),
    make_option("--mpi-host", default="", action="store", type="string",
                help="list of hosts to run job. Identical to --host for mpirun", metavar="HOSTS"),
    make_option("--mpi-host-file", default="", action="store", type="string",
                help="host list file to run job. Identical to --hostfile for mpirun", metavar="FILE"),
    make_option("--pbs", default=None, action="store_true",
                help="PBS for launching processes"),
    make_option("--pbs-name", default="cicada-pipe", action="store", type="string",
                help="PBS process name (up to 15 characters!) (default: cicada-pipe)", metavar="NAME"),
    make_option("--pbs-queue", default="ltg", action="store", type="string",
                help="PBS queue for launching processes (default: ltg)", metavar="NAME"),

    ## debug messages
    make_option("--debug", default=0, action="store", type="int"),
    ])

def run_command(command):
    try:
        retcode = subprocess.call(command, shell=True)
        if retcode:
            sys.exit(retcode)
    except:
        raise ValueError, "subprocess.call failed: %s" %(command)

def compressed_file(file):
    if not file:
        return file
    if os.path.exists(file):
        return file
    if os.path.exists(file+'.gz'):
	return file+'.gz'
    if os.path.exists(file+'.bz2'):
	return file+'.bz2'
    (base, ext) = os.path.splitext(file)
    if ext == '.gz' or ext == '.bz2':
	if os.path.exists(base):
	    return base
    return file

class Quoted:
    def __init__(self, arg):
        self.arg = arg
        
    def __str__(self):
        return '"' + str(self.arg) + '"'

class Option:
    def __init__(self, arg, value=None):
        self.arg = arg
        self.value = value

    def __str__(self,):
        option = self.arg
        
        if self.value is not None:
            if isinstance(self.value, int):
                option += " %d" %(self.value)
            elif isinstance(self.value, long):
                option += " %d" %(self.value)
            elif isinstance(self.value, float):
                option += " %.20g" %(self.value)
            else:
                option += " %s" %(str(self.value))
        return option

class Program:
    def __init__(self, *args):
        self.args = list(args[:])

    def __str__(self,):
        return ' '.join(map(str, self.args))
    
    def __iadd__(self, other):
        self.args.append(other)
        return self


class PBS:
    def __init__(self, queue=""):
        self.queue = queue
        self.qsub = 'qsub'
        
        # how to find binary location...?

    def run(self, command="", name="name", memory=0.0, mpi=None, threads=1, logfile=None):

        popen = subprocess.Popen("qsub -S /bin/sh", shell=True, stdin=subprocess.PIPE)
        pipe = popen.stdin
        
        pipe.write("#!/bin/sh\n")
        pipe.write("#PBS -S /bin/sh\n")
        pipe.write("#PBS -N %s\n" %(name))
        pipe.write("#PBS -W block=true\n")
        pipe.write("#PBS -e localhost:/dev/null\n")
        pipe.write("#PBS -o localhost:/dev/null\n")
        
        if self.queue:
            pipe.write("#PBS -q %s\n" %(self.queue))
            
        mem = ""
        if memory >= 1.0:
            mem=":mem=%dgb" %(int(memory))
        elif memory >= 0.001:
            mem=":mem=%dmb" %(int(memory * 1000))
        elif memory >= 0.000001:
            mem=":mem=%dkb" %(int(memory * 1000 * 1000))

        if mpi:
            pipe.write("#PBS -l select=%d:ncpus=%d:mpiprocs=1%s\n" %(mpi.number, threads, mem))
        else:
            pipe.write("#PBS -l select=1:ncpus=%d:mpiprocs=1%s\n" %(threads, mem))
        
        # setup variables
        if os.environ.has_key('TMPDIR_SPEC'):
            pipe.write("export TMPDIR_SPEC=%s\n" %(os.environ['TMPDIR_SPEC']))
        if os.environ.has_key('LD_LIBRARY_PATH'):
            pipe.write("export LD_LIBRARY_PATH=%s\n" %(os.environ['LD_LIBRARY_PATH']))
        if os.environ.has_key('DYLD_LIBRARY_PATH'):
            pipe.write("export DYLD_LIBRARY_PATH=%s\n" %(os.environ['DYLD_LIBRARY_PATH']))
        
        pipe.write("if test \"$PBS_O_WORKDIR\" != \"\"; then\n")
        pipe.write("  cd $PBS_O_WORKDIR\n")
        pipe.write("fi\n")

        prefix = ''
        if mpi:
            prefix = mpi.mpirun

            if os.environ.has_key('TMPDIR_SPEC'):
                prefix += ' -x TMPDIR_SPEC'
            if os.environ.has_key('LD_LIBRARY_PATH'):
                prefix += ' -x LD_LIBRARY_PATH'
            if os.environ.has_key('DYLD_LIBRARY_PATH'):
                prefix += ' -x DYLD_LIBRARY_PATH'
            prefix += ' '
        
        suffix = ''
        if logfile:
            suffix = " 2> %s" %(logfile)

        pipe.write(prefix + command + suffix + '\n')
        
        pipe.close()
        popen.wait()

class MPI:
    def __init__(self, dir="", hosts="", hosts_file="", number=0):
        
	self.dir = dir
	self.hosts = hosts
        self.hosts_file = hosts_file
        self.number = number
	
        if self.dir:
            if not os.path.exists(self.dir):
                raise ValueError, self.dir + " does not exist"
            self.dir = os.path.realpath(self.dir)

        if self.hosts_file:
            if not os.path.exists(self.hosts_file):
                raise ValueError, self.hosts_file + " does no exist"
            self.hosts_file = os.path.realpath(hosts_file)

        self.bindir = self.dir
	
        for binprog in ['mpirun']:
            if self.bindir:
                prog = os.path.join(self.bindir, 'bin', binprog)
                if not os.path.exists(prog) or os.path.isdir(prog):
                    prog = os.path.join(self.bindir, binprog)
                    if not os.path.exists(prog) or os.path.isdir(prog):
                        raise ValueError, prog + " does not exist at " + self.bindir
                    
                setattr(self, binprog, prog)
            else:
                setattr(self, binprog, binprog)
                
    def run(self, command, logfile=None):
        mpirun = self.mpirun
        if self.number > 0:
            mpirun += ' --np %d' %(self.number)
        if self.hosts:
            mpirun += ' --host %s' %(self.hosts)
        elif self.hosts_file:
            mpirun += ' --hostfile "%s"' %(self.hosts_file)

        if os.environ.has_key('TMPDIR_SPEC'):
            mpirun += ' -x TMPDIR_SPEC'
        if os.environ.has_key('LD_LIBRARY_PATH'):
            mpirun += ' -x LD_LIBRARY_PATH'
        if os.environ.has_key('DYLD_LIBRARY_PATH'):
            mpirun += ' -x DYLD_LIBRARY_PATH'

	mpirun += ' ' + command

        if logfile:
            mpirun += " 2> %s" %(logfile)
        
	run_command(mpirun)

class QSub:
    def __init__(self, mpi=None, pbs=None):
        self.mpi = mpi
        self.pbs = pbs
        
    def run(self, command, name="name", memory=0.0, threads=1, logfile=None):
        if self.pbs:
            self.pbs.run(str(command), name=name, memory=memory, threads=threads, logfile=logfile)
        else:
            if logfile:
                run_command(str(command) + " 2> %s" %(logfile))
            else:
                run_command(str(command))
    
    def mpirun(self, command, name="name", memory=0.0, threads=1, logfile=None):
        if not self.mpi:
            raise ValueError, "no mpi?"

        if self.pbs:
            self.pbs.run(str(command), name=name, memory=memory, mpi=self.mpi, logfile=logfile)
        else:
            self.mpi.run(str(command), logfile=logfile)

class CICADA:
    def __init__(self, dir=""):
        bindirs = []
        
        if not dir:
            dir = os.path.abspath(os.path.dirname(__file__))
            bindirs.append(dir)
            parent = os.path.dirname(dir)
            if parent:
                dir = parent
        else:
            dir = os.path.realpath(dir)
            if not os.path.exists(dir):
                raise ValueError, dir + " does not exist"
            bindirs.append(dir)
        
	for subdir in ('bin', 'progs', 'scripts'): 
	    bindir = os.path.join(dir, subdir)
	    if os.path.exists(bindir) and os.path.isdir(bindir):
		bindirs.append(bindir)
	
        for binprog in ('mpipe', ### mpi-launcher
                        'thrpe', ### thread-launcher
                        ):
	    
	    for bindir in bindirs:
		prog = os.path.join(bindir, binprog)
                
                if not os.path.exists(prog): continue
                if os.path.isdir(prog): continue
                
                setattr(self, binprog, prog)
                break

	    if not hasattr(self, binprog):
		raise ValueError, binprog + ' does not exist'

if __name__ == '__main__':
    (options, args) = opt_parser.parse_args()

    if not options.command:
        raise ValueError, "no command?"

    cicada = CICADA(dir=options.cicada_dir)
    
    ### MPI
    mpi = None
    if options.mpi_host or options.mpi_host_file or options.mpi > 0:
        mpi = MPI(dir=options.mpi_dir,
                  hosts=options.mpi_host,
                  hosts_file=options.mpi_host_file,
                  number=options.mpi)
    
    ### PBS
    pbs = None
    if options.pbs:
        pbs = PBS(queue=options.pbs_queue)
    
    ### QSUB
    qsub = QSub(mpi=mpi, pbs=pbs)

    if mpi:
        qsub.mpirun(Program(cicada.mpipe,
                            Option('--input', options.input),
                            Option('--output', options.output),
                            Option('--command', options.command)),
                    name=options.pbs_name,
                    memory=options.max_malloc,
                    threads=options.threads,
                    logfile=options.logfile)
    else:
        qsub.run(Program(cicada.thrpe,
                         Option('--input', options.input),
                         Option('--output', options.output),
                         Option('--command', options.command),
                         Option('--threads', options.threads)),
                 name=options.pbs_name,
                 memory=options.max_malloc,
                 threads=options.threads,
                 logfile=options.logfile)
