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
import shutil

import UserList
import UserString
import cStringIO

from optparse import OptionParser, make_option

opt_parser = OptionParser(
    option_list=[
    # output directory/filename prefix
    make_option("--root-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="root directory for outputs"),
    make_option("--corpus-dir", default="", action="store", type="string",
                metavar="PREFIX", help="corpus directory (default: ${root_dir}/corpus)"),
    make_option("--giza-f2e", default="", action="store", type="string",
                metavar="DIRECTORY", help="giza directory for P(f|e) (default: ${root_dir}/giza.${f}-${e})"),
    make_option("--giza-e2f", default="", action="store", type="string",
                metavar="DIRECTORY", help="giza directory for P(e|f) (default: ${root_dir}/giza.${e}-${f})"),
    make_option("--model-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="model directory (default: ${root_dir}/model)"),
    make_option("--alignment-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="alignment directory (default: ${model_dir})"),
    make_option("--lexical-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="lexical transltion table directory (default: ${model_dir)"),

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
    
    ## additional feature functions
    make_option("--lexicon-model1",             default=None, action="store_true", help="compute Model1 features"),
    make_option("--lexicon-noisy-or",           default=None, action="store_true", help="compute noisy-or features"),
    make_option("--lexicon-insertion-deletion", default=None, action="store_true", help="compute insertion/deletion features"),
    
    make_option("--threshold-insertion", default=0.01, action="store", type="float", help="threshold for insertion (default: 0.01)"),
    make_option("--threshold-deletion",  default=0.01, action="store", type="float", help="threshold for deletion (default: 0.01)"),
    
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

class QSUB(multiprocessing.Process):
    def __init__(self, command=""):
        
        multiprocessing.Process.__init__(self)
        self.command = command
        
    def run(self):
        popen = subprocess.Popen("qsub -S /bin/sh", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        popen.communicate(self.command)
        
class PBS:
    ### rewrite this by threads!
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
        for worker in self.workers:
            worker.join()
            
    def run(self, command="", threads=1, memory=0.0, name="name", logfile=None):
        pipe = cStringIO.StringIO()
        
        pipe.write("#!/bin/sh\n")
        pipe.write("#PBS -N %s\n" %(name))
        pipe.write("#PBS -e /dev/null\n")
        pipe.write("#PBS -o /dev/null\n")
        pipe.write("#PBS -W block=true\n")
        
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
        self.workers[-1].start();

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
    def __init__(self, cicada=None, model_dir="", cky=None):
        self.indexer = cicada.cicada_index_grammar
        self.filter  = cicada.cicada_filter_extract_phrase
        self.filter += " --cicada"
        self.cky = None
        self.grammar = "grammar"
        self.name = "phrase"
        
        self.counts = os.path.join(model_dir, "phrase-counts")
        self.scores = os.path.join(model_dir, "phrase-score")
        self.index  = os.path.join(model_dir, "phrase-index")
        
class IndexSCFG:
    def __init__(self, cicada=None, model_dir="", cky=None):
        self.indexer = cicada.cicada_index_grammar
        self.filter  = cicada.cicada_filter_extract_scfg
        self.filter += " --feature-root"
        self.cky = None
        self.grammar = "grammar"
        self.name = "scfg"

        self.counts = os.path.join(model_dir, "scfg-counts")
        self.scores = os.path.join(model_dir, "scfg-score")
        self.index  = os.path.join(model_dir, "scfg-index")

class IndexGHKM:
    def __init__(self, cicada=None, model_dir="", cky=None):
        self.indexer = cicada.cicada_index_tree_grammar
        self.filter  = cicada.cicada_filter_extract_ghkm
        self.cky = cky
        self.grammar = "tree-grammar"
        self.name = "ghkm"

        self.counts = os.path.join(model_dir, "ghkm-counts")
        self.scores = os.path.join(model_dir, "ghkm-score")
        self.index  = os.path.join(model_dir, "ghkm-index")


class IndexTree:
    def __init__(self, cicada=None, model_dir="", cky=None):
        self.indexer = cicada.cicada_index_tree_grammar
        self.filter  = cicada.cicada_filter_extract_ghkm
        self.cky = cky
        self.grammar = "tree-grammar"
        self.name = "tree"

        self.counts = os.path.join(model_dir, "tree-counts")
        self.scores = os.path.join(model_dir, "tree-score")
        self.index  = os.path.join(model_dir, "tree-index")

class Lexicon:
    def __init__(self, lexical_dir="",
                 model1=None,
                 noisy_or=None,
                 insertion_deletion=None,
                 threshold_insertion=0.01,
                 threshold_deletion=0.01):
        self.lexicon_source_target = compressed_file(os.path.join(lexical_dir, "lex.f2n"))
        self.lexicon_target_source = compressed_file(os.path.join(lexical_dir, "lex.n2f"))

        self.model1   = model1
        self.noisy_or = noisy_or
        self.insertion_deletion = insertion_deletion
        
        self.threshold_insertion = threshold_insertion
        self.threshold_deletion  = threshold_deletion


class Index(UserString.UserString):
    def __init__(self,
                 cicada=None,
                 indexer=None,
                 lexicon=None,
                 input="",
                 output="",
                 name="",
                 root_source="",
                 root_target="",
                 prior=0.1,
                 kbest=0,
                 quantize=None,
                 features=[],
                 attributes=[]):

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

        options_lexicon = ""
        if lexicon.model1 or lexicon.noisy_or or lexicon.insertion_deletion:
            options_lexicon += " --lexicon-source-target \"%s\"" %(lexicon.lexicon_source_target)
            options_lexicon += " --lexicon-target-source \"%s\"" %(lexicon.lexicon_target_source)
            
            if lexicon.model1:
                options_lexicon += " --model1"
            if lexicon.noisy_or:
                options_lexicon += " --noisy-or"
            if lexicon.insertion_deletion:
                options_lexicon += " --insertion-deletion"
                options_lexicon += " --threshold-insertion %.20g" %(lexicon.threshold_insertion)
                options_lexicon += " --threshold-deletion %.20g" %(lexicon.threshold_deletion)
        
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
            command += options_lexicon
            command += " | "
            command += indexer.indexer
        else:
            
            self.threads = 2

            command = indexer.filter
            command += " --dirichlet-prior %g" %(prior)
            command += " --root-source \"%s\"" %(root_source)
            command += " --root-target \"%s\"" %(root_target)
            command += options_lexicon
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
    def __init__(self, indexer=None):

        UserList.UserList.__init__(self)

        prefix = indexer.scores
        output = indexer.index

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

if options.root_dir:
    if not os.path.exists(options.root_dir):
	os.makedirs(options.root_dir)

if not options.model_dir:
    options.model_dir = os.path.join(options.root_dir, "model")
if not options.lexical_dir:
    options.lexical_dir = options.model_dir
if not options.alignment_dir:
    options.alignment_dir = options.model_dir

cicada = CICADA(options.cicada_dir)

indexer = None
if options.phrase:
    indexer = IndexPhrase(cicada, model_dir=options.model_dir, cky=options.cky)
elif options.scfg:
    indexer = IndexSCFG(cicada, model_dir=options.model_dir, cky=options.cky)
elif options.ghkm:
    indexer = IndexGHKM(cicada, model_dir=options.model_dir, cky=options.cky)
elif options.tree:
    indexer = IndexTree(cicada, model_dir=options.model_dir, cky=options.cky)
else:
    raise ValueError, "no indexer?"

scores = Scores(indexer)
lexicon = Lexicon(lexical_dir=options.lexical_dir,
                  model1=options.lexicon_model1,
                  noisy_or=options.lexicon_noisy_or,
                  insertion_deletion=options.lexicon_insertion_deletion,
                  threshold_insertion=options.threshold_insertion,
                  threshold_deletion=options.threshold_deletion)

fp = open(os.path.join(indexer.index, "files"), 'w')

if options.pbs:
    # we use pbs to run jobs
    pbs = PBS(queue=options.pbs_queue)
    
    for score in scores:
        index = Index(cicada=cicada,
                      indexer=indexer,
                      lexicon=lexicon,
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

        fp.write(os.path.basename(score.output)+'\n')

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
                      lexicon=lexicon,
                      input=score.input,
                      output=score.output,
                      root_source=scores.root_source,
                      root_target=scores.root_target,
                      prior=options.prior,
                      kbest=options.kbest,
                      quantize=options.quantize,
                      features=options.feature,
                      attributes=options.attribute)

        fp.write(os.path.basename(score.output)+'\n')

        mpi.run(command=index)
else:
    threads = Threads(cicada=cicada, threads=options.threads)
    
    for score in scores:
        index = Index(cicada=cicada,
                      indexer=indexer,
                      lexicon=lexicon,
                      input=score.input,
                      output=score.output,
                      root_source=scores.root_source,
                      root_target=scores.root_target,
                      prior=options.prior,
                      kbest=options.kbest,
                      quantize=options.quantize,
                      features=options.feature,
                      attributes=options.attribute)
        
        fp.write(os.path.basename(score.output)+'\n')

        threads.run(command=index)
