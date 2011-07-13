#!/usr/bin/env python
#
#  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#
###
### a wrapper script for indexing models, such as grammar, rule-grammar, word-clusters and global-lexicon
### Currently, we support only grammar and tree-grammar
###

import threading
import multiprocessing

import time
import sys
import os, os.path
import string
import re
import subprocess

import UserList
import UserString

from optparse import OptionParser, make_option

opt_parser = OptionParser(
    option_list=[
    make_option("--scores", default="", action="store", type="string", metavar="FILE", help="extracted scores"),
    make_option("--output", default="", action="store", type="string", metavar="FILE", help="output"),
    
    ## smoothing...
    make_option("--prior", default=0.1, action="store", type="float", metavar="PRIOR", help="model prior (default: 0.1)"),
    
    ## feature/attribute names
    make_option("--feature",   default=[], action="append", type="string", help="feature definitions"),
    make_option("--attribute", default=[], action="append", type="string", help="attribute definitions"),
    
    ### options...
    make_option("--phrase", default=None, action="store_true", help="index phrase grammar"),
    make_option("--scfg",   default=None, action="store_true", help="index SCFG grammar"),
    make_option("--ghkm",   default=None, action="store_true", help="index GHKM (tree-to-string) grammar"),
    make_option("--tree",   default=None, action="store_true", help="index tree-to-tree grammar"),
    
    make_option("--cky",    default=None, action="store_true", help="CKY mode indexing for tree-grammar"),
    
    ## quantize
    make_option("--quantize", default=None, action="store_true", help="perform quantization"),
    
    ## kbest options
    make_option("--kbest", default=0, action="store", type="float",
                metavar="KBEST", help="kbest max count of rules"),
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


### dump to stderr
stdout = sys.stdout
sys.stdout = sys.stderr

def run_command(command):
    fp = os.popen(command)
    while 1:
        data = fp.read(1)
        if not data: break
        stdout.write(data)

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
            
    def run(self, command="", threads=1, memory=0.0, name="name", logfile=None):
        popen = subprocess.Popen("qsub -S /bin/sh", shell=True, stdin=subprocess.PIPE)

        pipe = popen.stdin
        
        pipe.write("#!/bin/sh\n")
        pipe.write("#PBS -N %s\n" %(name))
        # we will run in non-block
        #pipe.write("#PBS -W block=true\n")
        pipe.write("#PBS -k n\n")
        
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
        
        pipe.close()
        popen.wait()

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
	
        for binprog in ('cicada_index_grammar',
                        'cicada_index_tree_grammar',
                        ### indexers
                        'cicada_filter_extract',
                        'cicada_filter_extract_phrase',
                        'cicada_filter_extract_scfg',
                        'cicada_filter_extract_ghkm',
                        ### filters
                        'mpish', ### mpi-launcher
                        'thrsh', ### thread-launcher
                        ### launchers
                        ):
	    
	    for bindir in self.bindirs:
		prog = os.path.join(bindir, binprog)
		if os.path.exists(prog):
		    setattr(self, binprog, prog)
		    break
	    if not hasattr(self, binprog):
		raise ValueError, binprog + ' does not exist'

class IndexPhrase:
    def __init__(self, cicada=None, cky=None):
        self.indexer = cicada.cicada_index_grammar
        self.filter  = cicada.cicada_filter_extract_phrase
        self.filter += " --cicada"
        self.cky = None
        self.name = "phrase"
        
class IndexSCFG:
    def __init__(self, cicada=None, cky=None):
        self.indexer = cicada.cicada_index_grammar
        self.filter  = cicada.cicada_filter_extract_scfg
        self.filter += " --feature-root"
        self.cky = None
        self.name = "scfg"

class IndexGHKM:
    def __init__(self, cicada=None, cky=None):
        self.indexer = cicada.cicada_index_grammar
        self.filter  = cicada.cicada_filter_extract_ghkm
        self.cky = cky
        self.name = "ghkm"

class IndexTree:
    def __init__(self, cicada=None, cky=None):
        self.indexer = cicada.cicada_index_grammar
        self.filter  = cicada.cicada_filter_extract_ghkm
        self.cky = cky
        self.name = "tree"

class Index(UserString.UserString):
    def __init__(self, cicada=None, indexer=None, input="", output="", name="", root_source="", root_target="", prior=0.1, kbest=0, quantize=None, features=[], attributes=[]):

        if not input:
            raise ValueError, "no input? %s" %(input)
        if not output:
            raise ValueError, "no output? %s" %(output)
        if not root_source:
            raise ValueError, "no root source? %s" %(root_source)
        if not root_target:
            raise ValueError, "no root target? %s" %(root_target)
        
        self.name    = "index-" + indexer.name
        self.logfile = "index-" + indexer.name + "." + name + ".log"
        
        command = ""

        if kbest > 0:
            
            self.threads = 3

            command = cicada.cicada_filter_extract
            command += " --nbest %d" %(kbest)
            command += " --input \"%s\"" %(input)
            command += " | "
            command += indexer.filter
            command += " --dirichlet-prior %g" %(prior)
            command += " --root-source \"%s\"" %(root_source)
            command += " --root-target \"%s\"" %(root_target)
            command += " | "
            command += indexer.indexer
        else:
            
            self.threads = 2

            command = indexer.filter
            command += " --dirichlet-prior %g" %(prior)
            command += " --root-source \"%s\"" %(root_source)
            command += " --root-target \"%s\"" %(root_target)
            command += " --input \"%s\"" %(input)
            command += " | "
            command += indexer.indexer
            
        if quantize:
            command += " --quantize"
        
        path = output
        sep = ':'
        if indexer.cky:
            path += sep
            sep = ','
            path += "cky=true"
        if features:
            path += sep
            sep = ','
            path += ','.join(features)
        if attributes:
            path += sep
            sep = ','
            path += ','.join(attributes)

        command += " --output \"%s\"" %(path)

        UserString.UserString.__init__(self, command)
        
class Score:
    def __init__(self, input="", output="", name=""):
        self.input = input
        self.output = output
        self.name = name


class Scores(UserList.UserList):
    def __init__(self, prefix="", output=""):

        UserList.UserList.__init__(self)

        path_files = os.path.join(prefix, 'files');
        path_root_source = os.path.join(prefix, 'root-source.gz')
        path_root_target = os.path.join(prefix, 'root-target.gz')
        
        if not os.path.exists(path_files):
            raise ValueError, "no path to files: %s" %(path_files)
        if not os.path.exists(path_root_source):
            raise ValueError, "no path to root-source: %s" %(path_root_source)
        if not os.path.exists(path_root_target):
            raise ValueError, "no path to root-target: %s" %(path_root_target)
        
        self.root_source = path_root_source
        self.root_target = path_root_target
        
        if output:
            if os.path.exists(output):
                if not os.path.isdir(output):
                    os.remove(output)
                    os.makedirs(output)
            else:
                os.makedirs(output)
        
        for line in open(path_files):
            name = line.strip()
            if not name: continue
            
            path = os.path.join(prefix, name)
            if not os.path.exists(path):
                raise ValueError, "no path to scores: %s" %(path)

            root,stem = os.path.splitext(name)

            self.append(Score(path, os.path.join(output, root + '.bin'), root))
            

(options, args) = opt_parser.parse_args()

cicada = CICADA(options.cicada_dir)

scores = Scores(options.scores, options.output)
indexer = None
if options.phrase:
    indexer = IndexPhrase(cicada, cky=options.cky)
elif options.scfg:
    indexer = IndexSCFG(cicada, cky=options.cky)
elif options.ghkm:
    indexer = IndexGHKM(cicada, cky=options.cky)
elif options.tree:
    indexer = IndexTree(cicada, cky=options.cky)
else:
    raise ValueError, "no indexer?"

if options.pbs:
    # we use pbs to run jobs
    pbs = PBS(queue=options.pbs_queue)
    
    for score in scores:
        index = Index(cicada=cicada,
                      indexer=indexer,
                      input=score.input,
                      output=score.output,
                      name=score.name,
                      root_source=scores.root_source,
                      root_target=scores.root_target,
                      prior=options.prior,
                      kbest=options.kbest,
                      quantize=options.quantize,
                      features=options.feature,
                      attributes=options.attribute)

        pbs.run(command=index, threads=index.threads, memory=options.max_malloc, name=index.name, logfile=index.logfile)
    
elif options.mpi:
    mpi = MPI(cicada=cicada,
              dir=options.mpi_dir,
              hosts=options.mpi_host,
              hosts_file=options.mpi_host_file,
              number=options.mpi)
        
    for score in scores:
        index = Index(cicada=cicada,
                      indexer=indexer,
                      input=score.input,
                      output=score.output,
                      root_source=scores.root_source,
                      root_target=scores.root_target,
                      prior=options.prior,
                      kbest=options.kbest,
                      quantize=options.quantize,
                      features=options.feature,
                      attributes=options.attribute)
        mpi.run(command=index)
else:
    threads = Threads(cicada=cicada, threads=options.threads)
    
    for score in scores:
        index = Index(cicada=cicada,
                      indexer=indexer,
                      input=score.input,
                      output=score.output,
                      root_source=scores.root_source,
                      root_target=scores.root_target,
                      prior=options.prior,
                      kbest=options.kbest,
                      quantize=options.quantize,
                      features=options.feature,
                      attributes=options.attribute)
        
        threads.run(command=index)
