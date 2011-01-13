#!/usr/bin/env python
#
#  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
#
### a wrapper script (similar to phrase-extract in moses)
### we support only "extraction" meaning only step 5 and 6
### TODO: use argparse for command-lines...?

import threading
import multiprocessing

import time
import sys
import os, os.path
import string
import re
import subprocess

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
    
    ### source/target flags
    make_option("--f", default="F", action="store", type="string",
                metavar="SUFFIX", help="source (or 'French')  language suffix for training corpus"),
    make_option("--e", default="E", action="store", type="string",
                metavar="SUFFIX", help="target (or 'English') language suffix for training corpus"),
    ### span...
    make_option("--sf", default="SF", action="store", type="string",
                metavar="SUFFIX", help="source (or 'French')  span suffix for training corpus"),
    make_option("--se", default="SE", action="store", type="string",
                metavar="SUFFIX", help="target (or 'English') span suffix for training corpus"),
    ### forest!
    make_option("--ff", default="SF", action="store", type="string",
                metavar="SUFFIX", help="source (or 'French')  forest suffix for training corpus"),
    make_option("--fe", default="SE", action="store", type="string",
                metavar="SUFFIX", help="target (or 'English') forest suffix for training corpus"),
    
    # data prefix
    make_option("--corpus", default="corpus", action="store", type="string",
                help="bilingual trainging corpus"),

    # alignment method
    make_option("--alignment", default="grow-diag-final-and", action="store", type="string",
                help="alignment methods (default: grow-diag-final-and)"),
    
    # steps
    make_option("--first-step", default=5, action="store", type="int", metavar='STEP', help="first step (default: 5)"),
    make_option("--last-step",  default=6, action="store", type="int", metavar='STEP', help="last step  (default: 6)"),


    ## option for lexicon
    make_option("--lexicon-prior", default=0.1, action="store", type="float", metavar="PRIOR", help="lexicon model prior (default: 0.1)"),
    
    # option for extraction
    make_option("--phrase", default=None, action="store_true", help="extract phrase"),
    make_option("--scfg",   default=None, action="store_true", help="extract SCFG"),
    make_option("--ghkm",   default=None, action="store_true", help="extract GHKM (tree-to-string)"),
    make_option("--tree",   default=None, action="store_true", help="extract tree-to-tree"),

    make_option("--non-terminal", default="[x]", action="store", type="string", help="default non-terminal for GHKM rule (default: [x])"),
    
    make_option("--max-span", default=15, action="store", type="int",
                metavar="LENGTH", help="maximum span size (default: 15)"),
    make_option("--min-hole", default=1, action="store", type="int",
                metavar="LENGTH", help="minimum hole size (default: 1)"),
    make_option("--max-length", default=7, action="store", type="int",
                metavar="LENGTH", help="maximum terminal length (default: 7)"),
    make_option("--max-fertility", default=4, action="store", type="int",
                metavar="FERTILITY", help="maximum terminal fertility (default: 4)"),
    make_option("--max-nodes", default=15, action="store", type="int",
                metavar="NODES", help="maximum rule nodes (default: 15)"),
    make_option("--max-height", default=4, action="store", type="int",
                metavar="HEIGHT", help="maximum rule height (default: 4)"),
    make_option("--exhaustive", default=None, action="store_true",
                help="exhaustive extraction in GHKM"),
    
    make_option("--ternary", default=None, action="store_true",
                help="extract ternary rule"),

    ## max-malloc
    make_option("--max-malloc", default=8, action="store", type="float",
                metavar="MALLOC", help="maximum memory in GB (default: 8)"),

    # CICADA Toolkit directory
    make_option("--toolkit-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="toolkit directory"),
    # MPI Implementation.. if different from standard location...
    make_option("--mpi-dir", default="", action="store", type="string",
                metavar="DIRECTORY", help="MPI directory"),

    # perform threading or MPI training    
    make_option("--mpi", default=0, action="store", type="int",
                help="# of processes for MPI-based parallel processing. Identical to --np for mpirun"),
    make_option("--mpi-host", default="", action="store", type="string",
                help="list of hosts to run job. Identical to --host for mpirun", metavar="HOSTS"),
    make_option("--mpi-host-file", default="", action="store", type="string",
                help="host list file to run job. Identical to --hostfile for mpirun", metavar="FILE"),
    
    make_option("--threads", default=2, action="store", type="int",
                help="# of thrads for thread-based parallel processing"),
    
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
            
    def run(self, command="", threads=1, memory=0.0, name="name", mpi=None, logfile=None):
        popen = subprocess.Popen("qsub -S /bin/sh", shell=True, stdin=subprocess.PIPE)

        pipe = popen.stdin
        
        pipe.write("#!/bin/sh\n")
        pipe.write("#PBS -N %s\n" %(name))
        pipe.write("#PBS -W block=true\n")
        
        if logfile:
            pipe.write("#PBS -e %s\n" %(logfile))
        else:
            pipe.write("#PBS -e /dev/null\n")
        pipe.write("#PBS -o /dev/null\n")
        
        if self.queue:
            pipe.write("#PBS -q %s\n" %(self.queue))

        if mpi:
            if memory > 0.0:
                if memory < 1.0:
                    pipe.write("#PBS -l select=%d:ncpus=3:mpiprocs=1:mem=%dmb\n" %(mpi.number, int(memory * 1000)))
                else:
                    pipe.write("#PBS -l select=%d:ncpus=3:mpiprocs=1:mem=%dgb\n" %(mpi.number, int(memory)))
            else:
                pipe.write("#PBS -l select=%d:ncpus=3:mpiprocs=1\n" %(mpi.number))
                
        else:
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

        if mpi:
            pipe.write("%s %s\n" %(mpi.mpirun, command))
        else:
            pipe.write("%s\n" %(command))
        
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
        if self.bindir:
            self.bindir = os.path.join(self.bindir, 'bin')
            if not os.path.exists(self.bindir) or not os.path.isdir(self.bindir):
                raise ValueError, self.bindir + " does not exist"
	
        for binprog in ['mpirun']:
            if self.bindir:
                prog = os.path.join(self.bindir, binprog)
                if not os.path.exists(prog):
                    raise ValueError, prog + " does not exist at " + self.bindir
                setattr(self, binprog, prog)
            else:
                setattr(self, binprog, binprog)
                
    def run(self, command):
        mpirun = self.mpirun
        #if self.dir:
        #    mpirun += ' --prefix %s' %(self.dir)
        if self.number > 0:
            mpirun += ' --np %d' %(self.number)
        if self.hosts:
            mpirun += ' --host %s' %(self.hosts)
        elif self.hosts_file:
            mpirun += ' --hostfile %s' %(self.hosts_file)
	mpirun += ' ' + command

	run_command(mpirun)


class Toolkit:
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

	if not self.bindirs:
	    raise ValueError, str(self.bindirs) + "  does not exist"
	
        for binprog in ('cicada_alignment',
                        ## step 4
                        'cicada_lexicon', 
                        ## step 5
                        'cicada_extract_phrase', 'cicada_extract_phrase_mpi',
                        'cicada_extract_scfg',   'cicada_extract_scfg_mpi',
                        'cicada_extract_ghkm',   'cicada_extract_ghkm_mpi',
                        'cicada_extract_tree',   'cicada_extract_tree_mpi',
                        ## step6
                        'cicada_extract_score', 'cicada_extract_score_mpi',):
	    
	    for bindir in self.bindirs:
		prog = os.path.join(bindir, binprog)
		if os.path.exists(prog):
		    setattr(self, binprog, prog)
		    break
	    if not hasattr(self, binprog):
		raise ValueError, binprog + ' does not exist'
        
class Corpus:

    def __init__(self, corpus="", f="", e="", sf="", se="", ff="", fe=""):

        self.source_tag = f
        self.target_tag = e

        self.source_span_tag = sf
        self.target_span_tag = se

        self.source_forest_tag = ff
        self.target_forest_tag = fe
        
        self.source = compressed_file(corpus+'.'+f)
        self.target = compressed_file(corpus+'.'+e)
        
        self.source_span = compressed_file(corpus+'.'+sf)
        self.target_span = compressed_file(corpus+'.'+se)
        
        self.source_forest = compressed_file(corpus+'.'+ff)
        self.target_forest = compressed_file(corpus+'.'+fe)

class Alignment:
    def __init__(self, alignment_dir="", alignment=""):
        self.alignment = compressed_file(os.path.join(alignment_dir, 'aligned.'+alignment))

        if not os.path.exists(self.alignment):
            raise ValueError, "no alignment data %s" %(self.alignment)
        

class Lexicon:
    def __init__(self, toolkit=None, corpus=None, alignment=None, lexical_dir="", prior=0.1,
                 threads=4, mpi=None, pbs=None,
                 debug=None):
        self.threads = threads
        self.mpi = mpi
        self.pbs = pbs
        
        self.source_target = compressed_file(os.path.join(lexical_dir, 'lex.f2n'))
        self.target_source = compressed_file(os.path.join(lexical_dir, 'lex.n2f'))
        self.makedirs = lexical_dir
        self.data = []

        self.data.append(corpus.source)
        self.data.append(corpus.target)

        command = "%s" %(toolkit.cicada_lexicon)
        
        command += " --source \"%s\"" %(corpus.source)
        command += " --target \"%s\"" %(corpus.target)
        command += " --alignment \"%s\"" %(alignment.alignment)
        
        command += " --output-source-target \"%s.gz\"" %(os.path.join(lexical_dir, 'lex.f2n'))
        command += " --output-target-source \"%s.gz\"" %(os.path.join(lexical_dir, 'lex.n2f'))
        
        command += " --variational-bayes"
        command += " --prior %g" %(prior)

        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

    def run(self):
        for file in self.data:
            if not os.path.exists(file):
                raise ValueError, "no file: " + file

        if not os.path.exists(self.makedirs):
            os.makedirs(self.makedirs)
        
        if self.pbs:
            self.pbs.run(command=self.command, threads=self.threads, name="lexicon", memory=8, logfile="lexicon.log")
        else:
            run_command(self.command)

        self.source_target = compressed_file(self.source_target)
        self.target_source = compressed_file(self.target_source)

class Extract:
    def __init__(self, max_malloc=8, threads=4, mpi=None, pbs=None, makedirs=""):
        self.threads = threads
        self.max_malloc = max_malloc
        self.mpi = mpi
        self.pbs = pbs
        self.makedirs = makedirs
        self.data = []
        self.logfile = ""
        self.name = ""
        
        if not hasattr(self, 'command'):
            self.command = ""
        
    def run(self):
        for file in self.data:
            if not os.path.exists(file):
                raise ValueError, "no file: " + file

        if not os.path.exists(self.makedirs):
            os.makedirs(self.makedirs)

        if not self.name:
            self.name = "extract"

        if self.mpi:
            if self.pbs:
                self.pbs.run(command=self.command, name=self.name, mpi=self.mpi, memory=self.max_malloc, logfile=self.logfile)
            else:
                self.mpi.run(self.command)
        else:
            if self.pbs:
                self.pbs.run(command=self.command, threads=self.threads, name="extract", memory=self.max_malloc, logfile=self.logfile)
            else:
                run_command(self.command)

class ExtractPhrase(Extract):
    
    def __init__(self, toolkit=None, corpus=None, alignment=None,
                 model_dir="",
                 max_length=7, max_fertility=4,
                 max_malloc=8, threads=4, mpi=None, pbs=None,
                 debug=None):
        Extract.__init__(self, max_malloc, threads, mpi, pbs, model_dir)
        
        self.counts = os.path.join(model_dir, "phrase-counts")

        self.data.append(corpus.source)
        self.data.append(corpus.target)
        self.logfile = "extract-phrase.log"
        self.name = "extract-phrase"
        
        prog_name = toolkit.cicada_extract_phrase
        if mpi:
            prog_name = toolkit.cicada_extract_phrase_mpi
        
        command = prog_name
        
        command += " --source \"%s\"" %(corpus.source)
        command += " --target \"%s\"" %(corpus.target)
        command += " --alignment \"%s\"" %(alignment.alignment)
        
        command += " --output \"%s\"" %(self.counts)
        
        command += " --max-length %d"    %(max_length)
        command += " --max-fertility %d" %(max_fertility)
        
        command += " --max-malloc %g" %(max_malloc)
        
        if not mpi:
            command += " --threads %d" %(threads)

        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

class ExtractSCFG(Extract):
    
    def __init__(self, toolkit=None, corpus=None, alignment=None,
                 model_dir="",
                 max_length=7, max_fertility=4, max_span=15, min_hole=1, ternary=None,
                 max_malloc=8, threads=4, mpi=None, pbs=None,
                 debug=None):
        Extract.__init__(self, max_malloc, threads, mpi, pbs, model_dir)
        
        self.counts = os.path.join(model_dir, "scfg-counts")
        
        self.data.append(corpus.source)
        self.data.append(corpus.target)
        self.logfile = "extract-scfg.log"
        self.name = "extract-scfg"
        
        if os.path.exists(corpus.source_span) and os.path.exists(corpus.target_span):
            raise ValueError, "both of source/target span specified... which one?"
        
        prog_name = toolkit.cicada_extract_scfg
        if mpi:
            prog_name = toolkit.cicada_extract_scfg_mpi
        
        command = prog_name
        
        command += " --source \"%s\"" %(corpus.source)
        command += " --target \"%s\"" %(corpus.target)
        command += " --alignment \"%s\"" %(alignment.alignment)

        if os.path.exists(corpus.source_span):
            command += " --spans-source \"%s\"" %(corpus.source_span)
        if os.path.exists(corpus.target_span):
            command += " --spans-target \"%s\"" %(corpus.target_span)
        
        command += " --output \"%s\"" %(self.counts)
        
        command += " --max-length %d"    %(max_length)
        command += " --max-fertility %d" %(max_fertility)
        command += " --max-span %d"      %(max_span)
        command += " --min-hole %d"      %(min_hole)

        if ternary:
            command += " --ternary"
        
        command += " --max-malloc %g" %(max_malloc)

        if not mpi:
            command += " --threads %d" %(self.threads)        
        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

class ExtractGHKM(Extract):
    
    def __init__(self, toolkit=None, corpus=None, alignment=None,
                 model_dir="",
                 non_terminal="", max_nodes=15, max_height=4,
                 exhaustive=None,
                 max_malloc=8, threads=4, mpi=None, pbs=None,
                 debug=None):
        Extract.__init__(self, max_malloc, threads, mpi, pbs, model_dir)
        
        self.counts = os.path.join(model_dir, "ghkm-counts")

        self.data.append(corpus.source_forest)
        self.data.append(corpus.target)
        self.logfile = "extract-ghkm.log"
        self.name = "extract-ghkm"
        
        prog_name = toolkit.cicada_extract_ghkm
        if mpi:
            prog_name = toolkit.cicada_extract_ghkm_mpi
        
        command = prog_name
        
        command += " --source \"%s\"" %(corpus.source_forest)
        command += " --target \"%s\"" %(corpus.target)
        command += " --alignment \"%s\"" %(alignment.alignment)
        
        command += " --output \"%s\"" %(self.counts)
        
        if non_terminal:
            if non_terminal[0] != '[' or non_terminal[-1] != ']':
                raise ValueError, "invalid non-terminal: %s" %(non_terminal)

            command += " --non-terminal \"%s\"" %(non_terminal)
        
        command += " --max-nodes %d"  %(max_nodes)
        command += " --max-height %d" %(max_height)

        if exhaustive:
            comamnd += " --exhaustive"
        
        command += " --max-malloc %g" %(max_malloc)

        if not mpi:
            command += " --threads %d" %(self.threads)
        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command


class ExtractTree(Extract):
    
    def __init__(self, toolkit=None, corpus=None, alignment=None,
                 model_dir="",
                 max_nodes=15, max_height=4,
                 max_malloc=8, threads=4, mpi=None, pbs=None,
                 debug=None):
        Extract.__init__(self, max_malloc, threads, mpi, pbs, model_dir)
        
        self.counts = os.path.join(model_dir, "tree-counts")

        self.data.append(corpus.source_forest)
        self.data.append(corpus.target_forest)
        self.logfile = "extract-tree.log"
        self.name = "extract-tree"
        
        prog_name = toolkit.cicada_extract_tree
        if mpi:
            prog_name = toolkit.cicada_extract_tree_mpi
        
        command = prog_name
        
        command += " --source \"%s\"" %(corpus.source_forest)
        command += " --target \"%s\"" %(corpus.target_forest)
        command += " --alignment \"%s\"" %(alignment.alignment)
        
        command += " --output \"%s\"" %(self.counts)
        
        command += " --max-nodes %d"  %(max_nodes)
        command += " --max-height %d" %(max_height)
        
        command += " --max-malloc %g" %(max_malloc)
        
        if not mpi:
            command += " --threads %d" %(self.threads)
            
        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"
        
        self.command = command

class ExtractScore(Extract):
    
    def __init__(self, toolkit=None, lexicon=None,
                 model_dir="",
                 phrase=None, scfg=None, ghkm=None, tree=None,
                 max_malloc=8, threads=4, mpi=None, pbs=None,
                 debug=None):
        Extract.__init__(self, max_malloc, threads, mpi, pbs, model_dir)
        
        option = ""
        if phrase:
            self.counts = os.path.join(model_dir, "phrase-counts")
            self.scores = os.path.join(model_dir, "phrase-score")
            option = " --score-phrase"
        elif scfg:
            self.counts = os.path.join(model_dir, "scfg-counts")
            self.scores = os.path.join(model_dir, "scfg-score")
            option = " --score-scfg"
        elif ghkm:
            self.counts = os.path.join(model_dir, "ghkm-counts")
            self.scores = os.path.join(model_dir, "ghkm-score")
            option = " --score-ghkm"
        elif tree:
            self.counts = os.path.join(model_dir, "tree-counts")
            self.scores = os.path.join(model_dir, "tree-score")
            option = " --score-ghkm"
        else:
            raise ValueError, "no count type?"

        if not os.path.exists(self.counts):
            raise ValueError, "no counts? %s" %(self.counts)

        self.logfile = "extract-score.log"
                
        prog_name = toolkit.cicada_extract_score
        if mpi:
            prog_name = toolkit.cicada_extract_score_mpi
        
        command = prog_name
        
        command += " --list \"%s\"" %(self.counts)
        command += " --output \"%s\"" %(self.scores)
        command += " --lexicon-source-target \"%s\"" %(lexicon.source_target)
        command += " --lexicon-target-source \"%s\"" %(lexicon.target_source)
        command += option
        command += " --max-malloc %g" %(max_malloc)
        
        if mpi:
            command += " --prog \"%s\"" %(prog_name)
        else:
            command += " --threads %d" %(self.threads)
            
        if debug:
            command += " --debug=%d" %(debug)
        else:
            command += " --debug"

        self.command = command


(options, args) = opt_parser.parse_args()

if not options.model_dir:
    options.model_dir = os.path.join(options.root_dir, "model")
if not options.lexical_dir:
    options.lexical_dir = options.model_dir
if not options.alignment_dir:
    options.alignment_dir = options.model_dir

toolkit = Toolkit(options.toolkit_dir)

mpi = None
if options.mpi_host or options.mpi_host_file or options.mpi > 0:
    mpi = MPI(dir=options.mpi_dir,
              hosts=options.mpi_host,
              hosts_file=options.mpi_host_file,
              number=options.mpi)

pbs = None
if options.pbs:
    pbs = PBS(queue=options.pbs_queue)

corpus = Corpus(corpus=options.corpus,
                f=options.f,
                e=options.e,
                sf=options.sf,
                se=options.se,
                ff=options.ff,
                fe=options.fe)

alignment = Alignment(options.alignment_dir, options.alignment)

lexicon = Lexicon(toolkit=toolkit, corpus=corpus, alignment=alignment,
                  lexical_dir=options.lexical_dir,
                  prior=options.lexicon_prior,
                  threads=options.threads, mpi=mpi, pbs=pbs,
                  debug=options.debug)

if options.first_step <= 4 and options.last_step >= 4:
    print "(4) generate lexical translation table started  @", time.ctime()
    lexicon.run()
    print "(4) generate lexical translation table finished @", time.ctime()

if options.first_step <= 5 and options.last_step >= 5:
    extract = None
    if options.phrase:
        extract = ExtractPhrase(toolkit=toolkit, corpus=corpus, alignment=alignment,
                                model_dir=options.model_dir,
                                max_length=options.max_length,
                                max_fertility=options.max_fertility,
                                max_malloc=options.max_malloc, threads=options.threads, mpi=mpi, pbs=pbs,
                                debug=options.debug)
    elif options.scfg:
        extract = ExtractSCFG(toolkit=toolkit, corpus=corpus, alignment=alignment,
                              model_dir=options.model_dir,
                              max_length=options.max_length,
                              max_fertility=options.max_fertility,
                              max_span=options.max_span,
                              min_hole=options.min_hole,
                              ternary=options.ternary,
                              max_malloc=options.max_malloc, threads=options.threads, mpi=mpi, pbs=pbs,
                              debug=options.debug)
    elif options.ghkm:
        extract = ExtractGHKM(toolkit=toolkit, corpus=corpus, alignment=alignment,
                              model_dir=options.model_dir,
                              non_terminal=options.non_terminal,
                              max_nodes=options.max_nodes,
                              max_height=options.max_height,
                              exhaustive=options.exhaustive,
                              max_malloc=options.max_malloc, threads=options.threads, mpi=mpi, pbs=pbs,
                              debug=options.debug)
    elif options.tree:
        extract = ExtractTree(toolkit=toolkit, corpus=corpus, alignment=alignment,
                              model_dir=options.model_dir,
                              max_nodes=options.max_nodes,
                              max_height=options.max_height,
                              max_malloc=options.max_malloc, threads=options.threads, mpi=mpi, pbs=pbs,
                              debug=options.debug)
    else:
        raise ValueError, "no count type?"

    print "(5) extract phrase table started @", time.ctime()
    extract.run()
    print "(5) extract phrase table finished @", time.ctime()

if options.first_step <= 6 and options.last_step >= 6:
    score = ExtractScore(toolkit=toolkit, lexicon=lexicon,
                         model_dir=options.model_dir,
                         phrase=options.phrase, scfg=options.scfg, ghkm=options.ghkm, tree=options.tree,
                         max_malloc=options.max_malloc, threads=options.threads, mpi=mpi, pbs=pbs,
                         debug=options.debug)
    
    print "(6) score phrase table started @", time.ctime()
    score.run()
    
    print "(6) score phrase table finished @", time.ctime()
