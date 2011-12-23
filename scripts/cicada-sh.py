#!/usr/bin/env python
#
#  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#
###
### a wrapper script for running multiple commands
### inspired by mpish and thrsh, we support pbs 
### Actually, we will run by mpish, thrsh and pbs
###

import threading
#import multiprocessing

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
    ## max-malloc
    make_option("--max-malloc", default=4, action="store", type="float",
                metavar="MALLOC", help="maximum memory in GB (default: 4)"),

    # CICADA Toolkit directory
    make_option("--cicada-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="cicada directory"),
    # MPI Implementation.. if different from standard location...
    make_option("--mpi-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="MPI directory"),

    make_option("--threads", default=2, action="store", type="int",
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
    make_option("--pbs-queue", default="ltg", action="store", type="string",
                help="PBS queue for launching processes (default: ltg)", metavar="NAME"),

    ## debug messages
    make_option("--debug", default=0, action="store", type="int"),
    ])

def run_command(command):
    retcode = subprocess.Popen(command, shell=True).wait()
    if retcode < 0:
        sys.exit(retcode)

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

class QSUB(threading.Process):
    def __init__(self, command=""):
        threading.Process.__init__(self)
        self.command = command
        self.qsub = None
        
    def run(self):
        popen = subprocess.Popen("qsub -S /bin/sh", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        data = popen.communicate(self.command)
        self.qsub = data[0].strip()
        
class PBS:
    def __init__(self, queue="", workingdir=os.getcwd()):
        self.queue = queue
        self.workingdir = workingdir
        self.tmpdir = None
        self.tmpdir_spec = None

        if os.environ.has_key('TMPDIR'):
            self.tmpdir = os.environ['TMPDIR']

        if os.environ.has_key('TMPDIR_SPEC'):
            self.tmpdir_spec = os.environ['TMPDIR_SPEC']

        self.workers = []

    def __del__(self):
        pass
        #for worker in self.workers:
        #    worker.join()
            
    def run(self, command="", threads=1, memory=0.0, name="cicada-sh", logfile=None):
        pipe = cStringIO.StringIO()
        
        pipe.write("#!/bin/sh\n")
        pipe.write("#PBS -N %s\n" %(name))
        pipe.write("#PBS -e /dev/null\n")
        pipe.write("#PBS -o /dev/null\n")
        #pipe.write("#PBS -W block=true\n")
        
        if self.workers and self.workers[-1].qsub:
            pipe.write("#PBS -W depend=after:%s\n" %(self.workers[-1].qsub))
        
        if self.queue:
            pipe.write("#PBS -q %s\n" %(self.queue))
            
        if memory > 0.0:
            if memory < 1.0:
                pipe.write("#PBS -l select=1:ncpus=%d:mpiprocs=1:mem=%dmb\n" %(threads, int(memory * 1000)))
            else:
                pipe.write("#PBS -l select=1:ncpus=%d:mpiprocs=1:mem=%dgb\n" %(threads, int(memory)))
        else:
            pipe.write("#PBS -l select=1:ncpus=%d:mpiprocs=1\n" %(threads))
        
        # setup TMPDIR and TMPDIR_SPEC
        if self.tmpdir:
            pipe.write("export TMPDIR=%s\n" %(self.tmpdir))
        if self.tmpdir_spec:
            pipe.write("export TMPDIR_SPEC=%s\n" %(self.tmpdir_spec))
            
        pipe.write("cd \"%s\"\n" %(self.workingdir))
        
        if logfile:
            pipe.write("%s >& %s\n" %(command, logfile))
        else:
            pipe.write("%s\n" %(command))

        self.workers.append(QSUB(pipe.getvalue()))
        self.workers[-1].start()
        self.workers[-1].join()

class Threads:
    
    def __init__(self, cicada=None, threads=1):
        command = "%s" %(cicada.thrsh)
        command += " --threads %d" %(threads)
        command += " --debug"
        
        self.popen = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE)
        self.pipe = self.popen.stdin
        
    def __del__(self):
        self.pipe.close()
        self.popen.wait()

    def run(self, command=""):
        self.pipe.write("%s\n" %(command))
        self.pipe.flush()

class MPI:
    
    def __init__(self, cicada=None, dir="", hosts="", hosts_file="", number=0):
        
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
                if not os.path.exists(prog):
                    prog = os.path.join(self.bindir, binprog)
                    if not os.path.exists(prog):
                        raise ValueError, prog + " does not exist at " + self.bindir
                    
                setattr(self, binprog, prog)
            else:
                setattr(self, binprog, binprog)
        
        command = self.mpirun
        if self.number > 0:
            command += ' --np %d' %(self.number)
            
        if self.hosts:
            command += ' --host %s' %(self.hosts)
        elif self.hosts_file:
            command += ' --hostfile %s' %(self.hosts_file)

        command += " %s" %(cicada.mpish)
        command += " --debug"
        
        self.popen = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE)
        self.pipe = self.popen.stdin
        
    def __del__(self):
        self.pipe.close()
        self.popen.wait()
        
    def run(self, command=""):
        self.pipe.write("%s\n" %(command))
        self.pipe.flush()


class CICADA:
    def __init__(self, dir=""):

	self.dir = dir	
	if not dir: return
	
	if not os.path.exists(self.dir):
	    raise ValueError, self.dir + " does not exist"
	
	self.dir = os.path.realpath(self.dir)
        
	self.bindirs = []
	for dir in ('bin', 'progs', 'scripts'): 
	    bindir = os.path.join(self.dir, dir)
	    if os.path.exists(bindir) and os.path.isdir(bindir):
		self.bindirs.append(bindir)
        self.bindirs.append(self.dir)
	
        for binprog in ('mpish', ### mpi-launcher
                        'thrsh', ### thread-launcher
                        ):
	    
	    for bindir in self.bindirs:
		prog = os.path.join(bindir, binprog)
		if os.path.exists(prog):
		    setattr(self, binprog, prog)
		    break
	    if not hasattr(self, binprog):
		raise ValueError, binprog + ' does not exist'


(options, args) = opt_parser.parse_args()

if options.pbs:
    # we use pbs to run jobs
    pbs = PBS(queue=options.pbs_queue)
    
    for line in sys.stdin:
        line = line.strip()
        if line:
            pbs.run(command=line, threads=options.threads, memory=options.max_malloc)

elif options.mpi:
    mpi = MPI(cicada=cicada,
              dir=options.mpi_dir,
              hosts=options.mpi_host,
              hosts_file=options.mpi_host_file,
              number=options.mpi)
    
    for line in sys.stdin:
        line = line.strip()
        if line:
            mpi.run(command=line)
else:
    threads = Threads(cicada=cicada, threads=options.threads)
    
    for line in sys.stdin:
        line = line.strip()
        if line:
            threads.run(command=line)
